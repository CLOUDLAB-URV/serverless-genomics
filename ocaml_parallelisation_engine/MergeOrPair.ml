(* This program is
    (C) 2019 Paolo Ribeca, <paolo.ribeca@gmail.com>.
   Please do not redistribute, or contact the author before doing so *)

(* This program only works for SE alignments - we should add a check *)

module Misc = struct
  (* Should be moved to Tools.ml *)

  let array_of_rlist l = Array.of_list (List.rev l)

  let array_riter f a =
    let l = Array.length a in
    let red_l = l - 1 in
    for i = red_l downto 0 do
      f a.(i)
    done

  let array_riteri f a =
    let l = Array.length a in
    let red_l = l - 1 in
    for i = red_l downto 0 do
      f i a.(i)
    done

  (* Should be moved to Tools.ml ? *)

  let input_line_num = ref 0

  let line__failwith s =
    Printf.eprintf "On line %d: %s\n%!" !input_line_num s;
    exit 1

  (* Should be moved to Tools.ml *)

  (* An ordered multimap is a map 'key -> 'val Set (no duplications allowed) *)
  module OrderedMultimap (OKey : Map.OrderedType) (OVal : Set.OrderedType) =
  struct
    (* Keys have type OKey.t, values OVal.t *)
    module KeyMap = Map.Make (OKey)
    module ValSet = Set.Make (OVal)

    (* To store the values *)
    type t = ValSet.t KeyMap.t

    let empty = KeyMap.empty
    let is_empty = KeyMap.is_empty

    let add k v om =
      try
        let s = KeyMap.find k om in
        KeyMap.add k (ValSet.add v s) (KeyMap.remove k om)
      with Not_found -> KeyMap.add k (ValSet.singleton v) om

    let remove_set = KeyMap.remove

    let remove k v om =
      try
        let s = KeyMap.find k om in
        let s = ValSet.remove v s in
        if ValSet.is_empty s then KeyMap.remove k om
        else KeyMap.add k (ValSet.remove v s) (KeyMap.remove k om)
      with Not_found -> om

    let iter_set = KeyMap.iter
    let iter f = KeyMap.iter (fun k s -> ValSet.iter (fun v -> f k v) s)
    let max_binding = KeyMap.max_binding
  end

  module TransitiveClosure = struct
    module IntSet = Set.Make (struct
      type t = int

      let compare = compare
    end)

    (* To make classes unique in the end *)
    module ClassSet = Set.Make (struct
      type t = IntSet.t

      let compare a b = if a == b then 0 else compare a b
    end)

    type t = IntSet.t array

    let create n =
      (* At the beginning all classes are different *)
      Array.init n (fun i -> IntSet.singleton i)

    let same_class tc i j =
      let l = Array.length tc in
      if i >= l || j >= l then
        failwith "TransitiveClosure.same_class: Invalid arguments";
      let merged = IntSet.union tc.(i) tc.(j) in
      IntSet.iter (fun id -> tc.(id) <- merged) merged

    let get_classes tc =
      let res = ref ClassSet.empty in
      Array.iter (fun cl -> res := ClassSet.add cl !res) tc;
      !res

    let iter f_begin_class f_elt f_end_class tc =
      ClassSet.iter
        (fun cl ->
          f_begin_class cl;
          IntSet.iter f_elt cl;
          f_end_class cl)
        (get_classes tc)
  end

  let re_summary = Str.regexp "[:+]"
  let re_separators = Str.regexp "[ \t]+"
  let re_alnum = Str.regexp "[0-9A-Za-z]"

  module Strand = struct
    type t = Forward | Reverse

    let parse failwith = function
      | "+" -> Forward
      | "-" -> Reverse
      | w ->
          failwith (Printf.sprintf "Expected strand, found '%s'" w);
          assert false (* Should never get here *)

    let unparse buf = function
      | Forward -> Buffer.add_char buf '+'
      | Reverse -> Buffer.add_char buf '-'

    let to_string str =
      let buf = Buffer.create 8 in
      unparse buf str;
      Buffer.contents buf
  end

  module IntMap = Map.Make (struct
    type t = int

    let compare = compare
  end)

  module IntSet = Set.Make (struct
    type t = int

    let compare = compare
  end)

  module StringMap = Map.Make (String)
end

(* This module provides a lightweight parser for the GEM format.
    It works on the record as an external string, and iterates over alignments implicitly
     by keeping track of the parsing state -- there could be millions of alignments for one read. *)
module GEMRecordParser = struct
  module Gigar = struct
    type op_t =
      | Match of int
      | Substitution of char
      | SkipQuery of int
      | SkipIndex of int * string
      | Splice of int * Misc.Strand.t option
      | Trim of int

    and t = op_t array

    let iter_intervals ?(mismatches_split_intervals = false) f_interv lo ops =
      let lo = ref lo and hi = ref (lo - 1) in
      Array.iter
        (function
          | Match l -> hi := !hi + l
          | w ->
              if
                (not mismatches_split_intervals)
                && match w with Substitution _ -> true | _ -> false
              then incr hi
              else (
                if !hi - !lo + 1 > 0 then f_interv !lo !hi;
                lo := !hi + 1;
                (match w with
                | Match _ -> assert false
                | Substitution _ -> incr lo
                | SkipIndex (l, _) | Splice (l, _) -> lo := !lo + l
                | SkipQuery _ | Trim _ -> ());
                hi := !lo - 1))
        ops;
      if !hi - !lo + 1 > 0 then f_interv !lo !hi

    let get_center_of_gravity lo ops =
      let num = ref 0 and den = ref 0 in
      iter_intervals
        (fun lo hi ->
          let len = hi - lo + 1 in
          num := !num + (len * (lo - 1)) + (len * (len + 1) / 2);
          den := !den + len)
        lo ops;
      int_of_float (float_of_int !num /. float_of_int !den)

    (* '>' ([+-]*[0-9]+) ([ACGNTacgnt]* ) ([-+*/%]) ) | '(' ([0-9]+) ')' | ([0-9]+) | ([ACGNTacgnt])
           ^             ^                ^                 ^              ^          ^
           1             2                3                 4              5          6 *)
    let gigar_re =
      Str.regexp
        ">\\([+-]*[0-9]+\\)\\([ACGNTacgnt]*\\)\\([-+*/%]\\)\\|(\\([0-9]+\\))\\|\\([0-9]+\\)\\|\\([ACGTNacgtn]\\)"

    let parse failwith strand s =
      let len = String.length s and idx = ref 0 and res = ref [] in
      while !idx < len do
        if not (Str.string_match gigar_re s !idx) then
          failwith (Printf.sprintf "Invalid element '%s' in GIGAR string" s);
        (*Printf.printf "FOUND, POS=%d, RES='%s'\n%!" !idx (Str.matched_string s);*)
        Tools.Misc.accum res
          (try
             (* We sort by the expected frequency *)
             Match (int_of_string (Str.matched_group 5 s))
           with Not_found -> (
             try
               match Str.matched_group 6 s with
               | "A" | "a" -> Substitution 'A'
               | "C" | "c" -> Substitution 'C'
               | "G" | "g" -> Substitution 'G'
               | "N" | "n" -> Substitution 'N'
               | "T" | "t" -> Substitution 'T'
               | _ -> assert false
             with Not_found -> (
               try
                 match Str.matched_group 3 s with
                 | "-" -> SkipQuery (int_of_string (Str.matched_group 1 s))
                 | "+" ->
                     let skip = int_of_string (Str.matched_group 1 s)
                     and seq =
                       try Str.matched_group 2 s with Not_found -> ""
                     in
                     let seq_len = String.length seq in
                     if skip = 0 then failwith "Zero skip in GIGAR string";
                     if seq_len > 0 && seq_len <> skip then
                       failwith "Unmatched skip and sequence in GIGAR string";
                     SkipIndex (skip, seq)
                 | "*" -> Splice (int_of_string (Str.matched_group 1 s), None)
                 | "/" ->
                     Splice (int_of_string (Str.matched_group 1 s), Some Forward)
                 | "%" ->
                     Splice (int_of_string (Str.matched_group 1 s), Some Reverse)
                 | _ -> assert false
               with Not_found -> (
                 try Trim (int_of_string (Str.matched_group 4 s))
                 with Not_found -> assert false))));
        idx := Str.match_end ()
      done;
      match strand with
      | Misc.Strand.Forward -> Misc.array_of_rlist !res
      | Reverse -> Array.of_list !res

    let unparse buf strand =
      (match strand with
      | Misc.Strand.Forward -> Array.iter
      | Reverse -> Misc.array_riter) (function
        | Match len -> Buffer.add_string buf (string_of_int len)
        | Substitution c -> Buffer.add_char buf c
        | SkipQuery len ->
            Buffer.add_char buf '>';
            Buffer.add_string buf (string_of_int len);
            Buffer.add_char buf '-'
        | SkipIndex (len, what) ->
            Buffer.add_char buf '>';
            Buffer.add_string buf (string_of_int len);
            Buffer.add_string buf what;
            Buffer.add_char buf '+'
        | Splice (len, dir) -> (
            Buffer.add_char buf '>';
            Buffer.add_string buf (string_of_int len);
            match dir with
            | None -> Buffer.add_char buf '*'
            | Some Misc.Strand.Forward -> Buffer.add_char buf '/'
            | Some Reverse -> Buffer.add_char buf '%')
        | Trim len ->
            Buffer.add_char buf '(';
            Buffer.add_string buf (string_of_int len);
            Buffer.add_char buf ')')

    let to_string strand ops =
      let buf = Buffer.create 256 in
      unparse buf strand ops;
      Buffer.contents buf

    type stats_t = {
      len_envelope : int;
      num_matches : int;
      len_matches : int;
      num_substitutions : int;
      num_query_skips : int;
      len_query_skips : int;
      num_index_skips : int;
      len_index_skips : int;
      num_splices : int;
      len_splices : int;
      num_trims : int;
      len_trims : int;
      num_errors : int;
      num_segments : int; (* Number of matches plus number of trims *)
    }

    let get_num_errors ops =
      let res = ref 0 in
      Array.iter
        (function
          | Match _ -> ()
          | Substitution _ | SkipQuery _ | SkipIndex _ | Splice _ | Trim _ ->
              incr res)
        ops;
      !res

    let get_stats ?(verbose = false) ops =
      if verbose then
        Printf.eprintf ">Getting stats for %s..."
          (to_string Misc.Strand.Forward ops);
      let res =
        ref
          {
            len_envelope = 0;
            num_matches = 0;
            len_matches = 0;
            num_substitutions = 0;
            num_query_skips = 0;
            len_query_skips = 0;
            num_index_skips = 0;
            len_index_skips = 0;
            num_splices = 0;
            len_splices = 0;
            num_trims = 0;
            len_trims = 0;
            num_errors = 0;
            num_segments = 0;
          }
      in
      Array.iter
        (function
          | Match l ->
              res :=
                {
                  !res with
                  num_matches = !res.num_matches + 1;
                  len_matches = !res.len_matches + l;
                }
          | Substitution _ ->
              res :=
                { !res with num_substitutions = !res.num_substitutions + 1 }
          | SkipQuery l ->
              res :=
                {
                  !res with
                  num_query_skips = !res.num_query_skips + 1;
                  len_query_skips = !res.len_query_skips + l;
                }
          | SkipIndex (l, _) ->
              res :=
                {
                  !res with
                  num_index_skips = !res.num_index_skips + 1;
                  len_index_skips = !res.len_index_skips + l;
                }
          | Splice (l, _) ->
              res :=
                {
                  !res with
                  num_splices = !res.num_splices + 1;
                  len_splices = !res.len_splices + l;
                }
          | Trim l ->
              res :=
                {
                  !res with
                  num_trims = !res.num_trims + 1;
                  len_trims = !res.len_trims + l;
                })
        ops;
      res :=
        {
          !res with
          len_envelope =
            !res.len_matches + !res.num_substitutions + !res.len_index_skips
            + !res.len_splices;
          num_errors =
            !res.num_substitutions + !res.num_query_skips + !res.num_index_skips
            + !res.num_splices + !res.num_trims;
          num_segments = !res.num_matches + !res.num_trims;
        };
      if verbose then
        Printf.eprintf
          " { envelope=%d, matches=(%d,%d), errors=%d, ratio=%.3g, \
           substitutions=%d, query skips=(%d,%d), index_skips=(%d,%d), \
           splices=(%d,%d), trims=(%d,%d) }\n\
           %!"
          !res.len_envelope !res.num_matches !res.len_matches !res.num_errors
          (float_of_int !res.len_matches /. float_of_int !res.len_envelope)
          !res.num_substitutions !res.num_query_skips !res.len_query_skips
          !res.num_index_skips !res.len_index_skips !res.num_splices
          !res.len_splices !res.num_trims !res.len_trims;
      !res

    (* Function to compare alignments by their stats. Better is bigger *)
    let compare_by_stats stats_1 stats_2 =
      let c = ref (compare stats_2.num_errors stats_1.num_errors) in
      (* +1 if 2 has more errors, i.e. 1 is better *)
      if !c = 0 then (
        c := compare stats_2.num_segments stats_1.num_segments;
        (* +1 if 2 has more segments, i.e. 1 is better *)
        if !c = 0 then (
          c := compare stats_1.len_matches stats_2.len_matches;
          (* +1 if 1 has more matched bases, i.e. 1 is better *)
          if !c = 0 then c := compare stats_2.len_envelope stats_1.len_envelope
            (* +1 if 2 has longer envelope, i.e. 1 is better *)));
      !c
    (*
          if stats_1.num_errors > stats_2.num_errors then (* 1 is worse *)
            -1
          else if stats_1.num_errors < stats_2.num_errors then
            1
          else
            let get_ratio stats = float_of_int stats.len_matches /. float_of_int stats.len_envelope in
            let ratio_1 = get_ratio stats_1 and ratio_2 = get_ratio stats_2 in
            if ratio_1 < ratio_2 then (* 1 is worse *)
              -1
            else if ratio_1 > ratio_2 then
              1
            else
              0
*)
  end

  module Alignment = struct
    type t = {
      reference : string;
      strand : Misc.Strand.t;
      position : int;
      (* There could be additional fields, but we do not care about them here *)
      errors : int;
      gigar : Gigar.t;
    }

    let parse failwith s =
      let fields = Tools.Split.as_array Tools.Split.colon_re s in
      let strand = Misc.Strand.parse failwith fields.(1) in
      let gigar = Gigar.parse failwith strand fields.(3) in
      {
        reference = fields.(0);
        strand;
        position = int_of_string fields.(2);
        errors = Gigar.get_num_errors gigar;
        gigar;
      }

    let unparse buf alignment =
      Buffer.add_string buf alignment.reference;
      Buffer.add_char buf ':';
      Misc.Strand.unparse buf alignment.strand;
      Buffer.add_char buf ':';
      Buffer.add_string buf (string_of_int alignment.position);
      Buffer.add_char buf ':';
      Gigar.unparse buf alignment.strand alignment.gigar

    let to_string alignment =
      let buf = Buffer.create 128 in
      unparse buf alignment;
      Buffer.contents buf
  end

  module Header = struct
    type t = {
      name : string;
      sequence : string;
      qualities : string;
      summary : int array;
    }

    let unparse buf alignment =
      Buffer.add_string buf alignment.name;
      Buffer.add_char buf '\t';
      Buffer.add_string buf alignment.sequence;
      Buffer.add_char buf '\t';
      Buffer.add_string buf alignment.qualities;
      Buffer.add_char buf '\t';
      Array.iteri
        (fun i errors ->
          if i > 0 then Buffer.add_char buf ':';
          Buffer.add_string buf (string_of_int errors))
        alignment.summary
  end

  type state_t = SBeginning | SAlignments of parsing_t | SEnd

  and parsing_t = {
    header : Header.t;
    error : int;
    remaining : int; (* Number of matches left with the given error *)
    pos : int; (* Position within string *)
  }

  and token_t = NHeader of Header.t | NAlignment of Alignment.t

  (* *)
  let parse_line_iter failwith line state =
    match !state with
    | SBeginning ->
        let len = String.length line
        and tabs = Array.make 4 0
        and tab = ref 0
        and pos = ref 0 in
        (try
           while !pos < len do
             if line.[!pos] = '\t' then (
               tabs.(!tab) <- !pos;
               incr tab;
               if !tab = 4 then raise Exit);
             incr pos
           done;
           failwith "Too few tabs in line"
         with Exit -> ());
        let header =
          let get_nth n =
            String.sub line (tabs.(n - 1) + 1) (tabs.(n) - tabs.(n - 1) - 1)
          in
          {
            Header.name = String.sub line 0 tabs.(0);
            sequence = get_nth 1;
            qualities = get_nth 2;
            summary =
              Array.map int_of_string
                (Tools.Split.as_array Misc.re_summary (get_nth 3));
          }
        in
        (state :=
           let next_pos = !pos + 1 in
           (* We have to cope with the case whereby the alignments are "-" *)
           if next_pos = len - 1 && line.[next_pos] = '-' then SEnd
           else
             SAlignments { header; error = -1; remaining = 0; pos = next_pos });
        NHeader header
    | SAlignments state_payload ->
        let len = String.length line and pos = ref state_payload.pos in
        (try
           while !pos < len do
             if line.[!pos] = ',' then raise Exit;
             incr pos
           done
         with Exit -> ());
        let summary = state_payload.header.summary in
        let max_error = Array.length summary - 1
        and error = ref state_payload.error
        and remaining = ref state_payload.remaining in
        while !remaining = 0 && !error < max_error do
          incr error;
          remaining := summary.(!error)
        done;
        if !error > max_error then
          failwith "Too many matches for the given summary";
        decr remaining;
        state :=
          if !pos = len then (* End of string *)
            SEnd
          else
            SAlignments
              {
                state_payload with
                error = !error;
                remaining = !remaining;
                pos = !pos + 1;
              };
        let alignment =
          Alignment.parse failwith
            (String.sub line state_payload.pos (!pos - state_payload.pos))
        in
        if !error <> alignment.errors then
          failwith "Inconsistency between summary and alignment";
        NAlignment alignment
    | SEnd -> raise Exit

  (* *)
  let parse failwith f_h f_a line =
    let state = ref SBeginning in
    let header =
      match parse_line_iter failwith line state with
      | NHeader header -> header
      | NAlignment _ -> assert false
    in
    f_h header;
    while !state <> SEnd do
      let alignment =
        match parse_line_iter failwith line state with
        | NHeader _ -> assert false
        | NAlignment alignment -> alignment
      in
      f_a header alignment
    done
end

(* A static class to hash names of reference sequences *)
module ReferenceHasher = struct
  let static_hash = ref Misc.StringMap.empty
  let inverse_static_hash = ref Misc.IntMap.empty

  let hash s =
    try Misc.StringMap.find s !static_hash
    with Not_found ->
      let idx = Misc.StringMap.cardinal !static_hash in
      static_hash := Misc.StringMap.add s idx !static_hash;
      inverse_static_hash := Misc.IntMap.add idx s !inverse_static_hash;
      idx

  let unhash idx = Misc.IntMap.find idx !inverse_static_hash

  (* We will hash the reference names as integers *)
  type stranded_hashed_reference_t = {
    reference : hashed_reference_t;
    strand : Misc.Strand.t;
  }

  and hashed_reference_t = int

  (* A map holding the segments/intervals from the alignments.
     We will also hash the alignments as integers *)
  module StrandedHashedReferenceMap = Map.Make (struct
    type t = stranded_hashed_reference_t

    let compare = compare
  end)
end

(* Merges alignments from different records by detecting their overlaps *)
module AlignmentMerger = struct
  type interval_t =
    | BeginSegment of hashed_alignment_t
    | EndSegment of hashed_alignment_t

  and hashed_alignment_t = int

  (* To store the relation position on the stranded reference -> one or more interval extrema *)
  module IntOrderedMultimap =
    Misc.OrderedMultimap
      (struct
        type t = int

        let compare = compare
      end)
      (struct
        type t = interval_t

        let compare = compare
      end)

  type overlap_detector_t =
    IntOrderedMultimap.t ReferenceHasher.StrandedHashedReferenceMap.t

  (* A module to sort alignments and automatically get the best one *)
  module StatsMap = Map.Make (struct
    type t = GEMRecordParser.Gigar.stats_t

    let compare = GEMRecordParser.Gigar.compare_by_stats
  end)

  (* A module to store several alignments sorted by error -- will hold the result.
     Note that by now alignments must have already been made unique, so we can store them in an OrderedMultimap *)
  module SingleEndAlignmentsByError =
    Misc.OrderedMultimap
      (struct
        type t = int

        let compare = compare
      end)
      (struct
        type t = GEMRecordParser.Alignment.t

        let compare = compare
      end)

  let merge ?(verbose = false) failwith buf lines =
    let header = ref None in
    let hashed_alignments = ref Misc.IntMap.empty in
    let add_to_hashed_alignments alignment =
      let idx = Misc.IntMap.cardinal !hashed_alignments in
      hashed_alignments :=
        Misc.IntMap.add idx
          ( alignment,
            GEMRecordParser.Gigar.get_stats ~verbose
              alignment.GEMRecordParser.Alignment.gigar )
          !hashed_alignments;
      idx
    in
    let overlap_detector =
      ref ReferenceHasher.StrandedHashedReferenceMap.empty
    in
    let add_to_overlap_detector alignment =
      let hashed_alignment = add_to_hashed_alignments alignment in
      let hashed_reference = ReferenceHasher.hash alignment.reference in
      let stranded_hashed_reference =
        {
          ReferenceHasher.reference = hashed_reference;
          strand = alignment.strand;
        }
      in
      let interval_multimap =
        ref
          (try
             ReferenceHasher.StrandedHashedReferenceMap.find
               stranded_hashed_reference !overlap_detector
           with Not_found ->
             let interval_multimap = IntOrderedMultimap.empty in
             overlap_detector :=
               ReferenceHasher.StrandedHashedReferenceMap.add
                 stranded_hashed_reference interval_multimap !overlap_detector;
             interval_multimap)
      in
      (* We iterate over all the segment in the alignment, and add them to the overlap detector *)
      GEMRecordParser.Gigar.iter_intervals ~mismatches_split_intervals:true
        (fun lo_seg hi_seg ->
          if verbose then
            Printf.eprintf
              ">Adding segment (%d-%d) for alignment #%d, %s...\n%!" lo_seg
              hi_seg hashed_alignment
              (GEMRecordParser.Alignment.to_string alignment);
          interval_multimap :=
            IntOrderedMultimap.add lo_seg (BeginSegment hashed_alignment)
              !interval_multimap;
          interval_multimap :=
            IntOrderedMultimap.add hi_seg (EndSegment hashed_alignment)
              !interval_multimap)
        alignment.position alignment.gigar;
      (* We replace the interval map *)
      overlap_detector :=
        ReferenceHasher.StrandedHashedReferenceMap.add stranded_hashed_reference
          !interval_multimap !overlap_detector
    in
    Array.iter
      (fun line ->
        GEMRecordParser.parse Misc.line__failwith
          (fun h ->
            (* Apart from the summary, all headers must be the same *)
            match !header with
            | None -> header := Some { h with summary = [||] }
            | Some hh ->
                if { h with summary = [||] } <> hh then
                  failwith "Inconsistent alignment headers")
          (fun header -> add_to_overlap_detector)
          line)
      lines;
    (* We extract equivalence classes *)
    let classes =
      Misc.TransitiveClosure.create (Misc.IntMap.cardinal !hashed_alignments)
    in
    ReferenceHasher.StrandedHashedReferenceMap.iter
      (fun { reference; strand } interval_multimap ->
        (* Each stranded reference is processed separately *)
        let open_als = ref Misc.IntSet.empty in
        IntOrderedMultimap.iter_set
          (fun pos multimap ->
            (* We first process the Begin-s, and then the End-s (inclusive, as per UCSC's convention) *)
            IntOrderedMultimap.ValSet.iter
              (function
                | BeginSegment hashed_alignment ->
                    if verbose && !open_als <> Misc.IntSet.empty then
                      Printf.eprintf
                        ">@%s:%s:%d: Alignment #%d is in the same class as \
                         #%d...\n\
                         %!"
                        (ReferenceHasher.unhash reference)
                        (Misc.Strand.to_string strand)
                        pos hashed_alignment
                        (Misc.IntSet.min_elt !open_als);
                    if !open_als <> Misc.IntSet.empty then
                      Misc.TransitiveClosure.same_class classes
                        (Misc.IntSet.min_elt !open_als)
                        hashed_alignment;
                    open_als := Misc.IntSet.add hashed_alignment !open_als
                | EndSegment _ -> ())
              multimap;
            IntOrderedMultimap.ValSet.iter
              (function
                | BeginSegment _ -> ()
                | EndSegment hashed_alignment ->
                    open_als := Misc.IntSet.remove hashed_alignment !open_als)
              multimap)
          interval_multimap;
        assert (!open_als = Misc.IntSet.empty))
      !overlap_detector;
    (* We extract alignments from classes.
       Alignments are sorted by error *)
    let res = ref SingleEndAlignmentsByError.empty
    and current_class = ref StatsMap.empty in
    Misc.TransitiveClosure.iter
      (fun _ ->
        if verbose then Printf.eprintf ">Processing class...\n%!";
        current_class := StatsMap.empty)
      (fun hashed_alignment ->
        let alignment, stats =
          Misc.IntMap.find hashed_alignment !hashed_alignments
        in
        if verbose then
          Printf.eprintf ">>Adding alignment %s...\n%!"
            (GEMRecordParser.Alignment.to_string alignment);
        current_class := StatsMap.add stats hashed_alignment !current_class)
      (fun _ ->
        let _, hashed_best_alignment = StatsMap.max_binding !current_class in
        let best_alignment, _ =
          Misc.IntMap.find hashed_best_alignment !hashed_alignments
        in
        if verbose then
          Printf.eprintf ">>Best representative is %s.\n%!"
            (GEMRecordParser.Alignment.to_string best_alignment);
        res :=
          SingleEndAlignmentsByError.add best_alignment.errors best_alignment
            !res)
      classes;
    (* We de-parse alignments *)
    let buf_alignments = Buffer.create 4096 in
    (* First we reconstruct the alignment summary *)
    let max_errors =
      try
        let max_errors, _ = SingleEndAlignmentsByError.max_binding !res in
        max_errors
      with Not_found -> 0
    in
    let summary = Array.make (max_errors + 1) 0 in
    SingleEndAlignmentsByError.iter_set
      (fun errors alignments ->
        summary.(errors) <-
          SingleEndAlignmentsByError.ValSet.cardinal alignments;
        SingleEndAlignmentsByError.ValSet.iter
          (fun alignment ->
            if Buffer.length buf_alignments > 0 then
              Buffer.add_char buf_alignments ',';
            GEMRecordParser.Alignment.unparse buf_alignments alignment)
          alignments)
      !res;
    GEMRecordParser.Header.unparse buf
      { (Tools.Misc.unbox_opt !header) with GEMRecordParser.Header.summary };
    Buffer.add_char buf '\t';
    if Buffer.length buf_alignments > 0 then
      Buffer.add_buffer buf buf_alignments
    else Buffer.add_char buf '-';
    Buffer.add_char buf '\n'
end

(* Pairs alignments from two different records, according to distance and orientation *)
module AlignmentPairer = struct
  module PairedOrientations = struct
    type t = FF | FR | RF | RR

    let of_string = function
      | "FF" | "ff" | "++" -> FF
      | "FR" | "fr" | "+-" -> FR
      | "RF" | "rf" | "-+" -> RF
      | "RR" | "rr" | "--" -> RR
      | w ->
          Printf.sprintf
            "PairedOrientations.of_string: Unknown orientation '%s'" w
          |> failwith

    let to_string = function FF -> "FF" | FR -> "FR" | RF -> "RF" | RR -> "RR"
  end

  module UnpairedPrintMode = struct
    type t = Discard | Joint | Separate

    let of_string = function
      | "0" -> Discard
      | "1" -> Joint
      | "2" -> Separate
      | w ->
          Printf.sprintf "PairedPrintMode.of_string: Unknown print mode '%s'" w
          |> failwith

    let to_string = function Discard -> "0" | Joint -> "1" | Separate -> "2"
  end

  let find_common_name ?(accept_unpairable_names = false) failwith name_1 name_2
      =
    let error () =
      if accept_unpairable_names then name_1
      else (
        Printf.sprintf
          "AlignmentPairer.find_common_name: Unable to deduce common name for \
           the pair (names='%s','%s')"
          name_1 name_2
        |> failwith;
        assert false (* Not supposed to get here in the first place *))
    in
    (* First we split at spaces *)
    let split_at_separators s =
      (Tools.Split.as_array Misc.re_separators s).(0)
    in
    let name_1 = split_at_separators name_1
    and name_2 = split_at_separators name_2 in
    let l = String.length name_1 in
    if l <> String.length name_2 || l < 3 then error ()
    else
      let i = ref 0 in
      (try
         while !i < l do
           if name_1.[!i] = name_2.[!i] then incr i else raise Exit
         done
       with Exit -> ());
      let red_l = l - 1 in
      if !i >= l - 1 then
        if !i = l then (* Same names *)
          name_1
        else if name_1.[red_l] = '1' && name_2.[red_l] = '2' then
          if !i = red_l then String.sub name_1 0 red_l
          else
            (* Names such as PREFIX '/1' *)
            let red_red_l = red_l - 1 in
            if Str.string_match Misc.re_alnum (String.sub name_1 red_red_l 1) 0
            then String.sub name_1 0 red_red_l
            else error ()
        else error ()
      else error ()

  (* A module to store several alignments sorted by error -- will hold the result.
     Note that by now alignments must have already been made unique, so we can store them in an OrderedMultimap *)
  module PairedEndAlignmentsByError =
    Misc.OrderedMultimap
      (struct
        type t = int

        let compare = compare
      end)
      (struct
        type t = GEMRecordParser.Alignment.t * GEMRecordParser.Alignment.t

        let compare = compare
      end)

  let pair ?(verbose = false) ?(accept_unpairable_names = false)
      ?(unpaired_print_mode = UnpairedPrintMode.Joint) max_pairing_distance
      allowed_pairing_orientations max_pairing_number failwith buf lines =
    (* There must be exactly two inputs *)
    if Array.length lines <> 2 then
      Printf.sprintf "AlignmentPairer.merge: 2 inputs required, found %d"
        (Array.length lines)
      |> failwith;
    if verbose then (
      Printf.eprintf ">[1]: %s\n%!" lines.(0);
      Printf.eprintf ">[2]: %s\n%!" lines.(1));
    (* Alas, we have to fully parse the second line, otherwise we would spend a very long time parsing it over and over again *)
    let header = ref None and alignments_2 = ref [] in
    GEMRecordParser.parse failwith
      (fun header_2 -> header := Some header_2)
      (fun _ -> Tools.Misc.accum alignments_2)
      lines.(1);
    let headers_checked = ref false
    and num_paired = ref 0
    and res = ref PairedEndAlignmentsByError.empty in
    try
      GEMRecordParser.parse failwith
        (fun header_1 ->
          if not !headers_checked then (
            let header_2 = Tools.Misc.unbox_opt !header in
            header :=
              Some
                {
                  name =
                    find_common_name ~accept_unpairable_names failwith
                      header_1.name header_2.name;
                  sequence =
                    Printf.sprintf "%s %s" header_1.sequence header_2.sequence;
                  qualities =
                    Printf.sprintf "%s %s" header_1.qualities header_2.qualities;
                  summary = [||];
                };
            headers_checked := true))
        (fun _ alignment_1 ->
          let center_of_gravity_1 =
            GEMRecordParser.Gigar.get_center_of_gravity alignment_1.position
              alignment_1.gigar
          in
          List.iter
            (fun alignment_2 ->
              if verbose then
                Printf.eprintf ">>Trying to pair alignments %s and %s...%!"
                  (GEMRecordParser.Alignment.to_string alignment_1)
                  (GEMRecordParser.Alignment.to_string alignment_2);
              if
                (* References must be the same *)
                alignment_1.reference = alignment_2.reference
                (* Strands must be compatible *)
                && (match (alignment_1.strand, alignment_2.strand) with
                   | Misc.Strand.Forward, Misc.Strand.Forward ->
                       allowed_pairing_orientations.(0)
                   | Misc.Strand.Forward, Misc.Strand.Reverse ->
                       allowed_pairing_orientations.(1)
                   | Misc.Strand.Reverse, Misc.Strand.Forward ->
                       allowed_pairing_orientations.(2)
                   | Misc.Strand.Reverse, Misc.Strand.Reverse ->
                       allowed_pairing_orientations.(3))
                &&
                (* As for positions, we compare centers-of-gravity *)
                let center_of_gravity_2 =
                  GEMRecordParser.Gigar.get_center_of_gravity
                    alignment_2.position alignment_2.gigar
                in
                abs (center_of_gravity_1 - center_of_gravity_2)
                <= max_pairing_distance
              then (
                if verbose then Printf.eprintf " OK.\n%!";
                incr num_paired;
                if !num_paired <= max_pairing_number then
                  res :=
                    PairedEndAlignmentsByError.add
                      (alignment_1.errors + alignment_2.errors)
                      (alignment_1, alignment_2) !res)
              else if verbose then
                let center_of_gravity_2 =
                  GEMRecordParser.Gigar.get_center_of_gravity
                    alignment_2.position alignment_2.gigar
                in
                Printf.eprintf
                  " KO (%b %b %b, dist=%d, cogs=%d,%d, mpd=%d).\n%!"
                  (alignment_1.reference = alignment_2.reference)
                  (match (alignment_1.strand, alignment_2.strand) with
                  | Misc.Strand.Forward, Misc.Strand.Forward ->
                      allowed_pairing_orientations.(0)
                  | Misc.Strand.Forward, Misc.Strand.Reverse ->
                      allowed_pairing_orientations.(1)
                  | Misc.Strand.Reverse, Misc.Strand.Forward ->
                      allowed_pairing_orientations.(2)
                  | Misc.Strand.Reverse, Misc.Strand.Reverse ->
                      allowed_pairing_orientations.(3))
                  (abs (center_of_gravity_1 - center_of_gravity_2)
                  <= max_pairing_distance)
                  (abs (center_of_gravity_1 - center_of_gravity_2))
                  center_of_gravity_1 center_of_gravity_2 max_pairing_number)
            !alignments_2)
        lines.(0);
      (* Everything went fine, and we can output the paired records *)
      (* We de-parse alignments *)
      let buf_alignments = Buffer.create 1048576 in
      (* First we reconstruct the alignment summary *)
      let max_errors =
        try
          let max_errors, _ = PairedEndAlignmentsByError.max_binding !res in
          max_errors
        with Not_found -> 0
      in
      let summary = Array.make (max_errors + 1) 0 in
      PairedEndAlignmentsByError.iter_set
        (fun errors alignments ->
          summary.(errors) <-
            PairedEndAlignmentsByError.ValSet.cardinal alignments;
          PairedEndAlignmentsByError.ValSet.iter
            (fun (alignment_1, alignment_2) ->
              if Buffer.length buf_alignments > 0 then
                Buffer.add_char buf_alignments ',';
              GEMRecordParser.Alignment.unparse buf_alignments alignment_1;
              Buffer.add_string buf_alignments "::";
              GEMRecordParser.Alignment.unparse buf_alignments alignment_2)
            alignments)
        !res;
      if Buffer.length buf_alignments > 0 then (
        (* We always print paired alignments *)
        GEMRecordParser.Header.unparse buf
          { (Tools.Misc.unbox_opt !header) with GEMRecordParser.Header.summary };
        Buffer.add_char buf '\t';
        Buffer.add_buffer buf buf_alignments;
        Buffer.add_char buf '\n')
      else (
        if verbose then Printf.eprintf ">>No paired alignments.\n%!";
        match unpaired_print_mode with
        | UnpairedPrintMode.Discard -> ()
        | Joint ->
            GEMRecordParser.Header.unparse buf
              {
                (Tools.Misc.unbox_opt !header) with
                GEMRecordParser.Header.summary = [| 0 |];
              };
            Buffer.add_string buf "\t-\n"
        | Separate ->
            Buffer.add_string buf lines.(0);
            Buffer.add_char buf '\n';
            Buffer.add_string buf lines.(1);
            Buffer.add_char buf '\n')
    with Exit -> (
      (* In this case the number of paired alignments exceeded user-selected parameters,
          and we just bail out printing unpaired reads with changed names in order to signal what happened *)
      if verbose then Printf.eprintf ">>Too many paired alignments.\n%!";
      match unpaired_print_mode with
      | UnpairedPrintMode.Discard -> ()
      | Joint ->
          let header = Tools.Misc.unbox_opt !header in
          GEMRecordParser.Header.unparse buf
            {
              header with
              GEMRecordParser.Header.summary = [| 0 |];
              name = "!" ^ header.name;
            };
          Buffer.add_string buf "\t-\n"
      | Separate ->
          Buffer.add_char buf '!';
          Buffer.add_string buf lines.(0);
          Buffer.add_char buf '\n';
          Buffer.add_char buf '!';
          Buffer.add_string buf lines.(1);
          Buffer.add_char buf '\n')
end

module Defaults = struct
  let lines_per_block = 10000

  (* Only for pairing *)
  let max_pairing_distance = 500000
  let allowed_pairing_orientations = [| false; true; true; false |]
  let max_pairing_number = 10000
  let accept_unpairable_names = false
  let unpaired_print_mode = AlignmentPairer.UnpairedPrintMode.Joint

  (* *)
  let chunk_size = 64
  let threads = 1
  let verbose = false
end

module Parameters = struct
  module WorkingMode = struct
    type t = Merge | Pair

    let of_string = function
      | "merge" -> Merge
      | "pair" -> Pair
      | _ -> assert false
  end

  let working_mode = ref (None : WorkingMode.t option)
  let lines_per_block = ref Defaults.lines_per_block
  let max_pairing_distance = ref Defaults.max_pairing_distance
  let allowed_pairing_orientations = Defaults.allowed_pairing_orientations
  let max_pairing_number = ref Defaults.max_pairing_number
  let accept_unpairable_names = ref Defaults.accept_unpairable_names
  let unpaired_print_mode = ref Defaults.unpaired_print_mode
  let inputs = ref []
  let output = ref ""
  let chunk_size = ref Defaults.chunk_size
  let threads = ref Defaults.threads
  let verbose = ref Defaults.verbose
end

let version = "0.1"

let () =
  Printf.eprintf "This is the MergeOrPair program (version %s)\n%!" version;
  Printf.eprintf " (c) 2019 Paolo Ribeca, <paolo.ribeca@gmail.com>\n%!";
  let module TA = Tools.Argv in
  TA.parse
    [
      ([], None, [ "Working mode" ], TA.Optional, fun _ -> ());
      ( [ "merge"; "pair" ],
        None,
        [ "run in merge or pair mode" ],
        TA.Mandatory,
        fun mode ->
          Parameters.working_mode :=
            Some (Parameters.WorkingMode.of_string mode) );
      ([], None, [ "Pairing" ], TA.Optional, fun _ -> ());
      ( [ "-d"; "--maximum-pairing-distance" ],
        Some "<non_negative_integer>",
        [ "maximum index distance for a pairing to be made" ],
        TA.Default (fun () -> string_of_int Defaults.max_pairing_distance),
        fun _ ->
          Parameters.max_pairing_distance := TA.get_parameter_int_non_neg () );
      ( [ "--set-orientation"; "--unset-orientation" ],
        Some "FF|FR|RF|RR",
        [ "allow/disallow pairings with the given orientation" ],
        TA.Default
          (fun () ->
            let buf = Buffer.create 128 in
            Buffer.add_char buf '{';
            Array.iteri
              (fun i b ->
                Buffer.add_string buf
                  (match i with
                  | 0 -> " FF="
                  | 1 -> " FR="
                  | 2 -> " RF="
                  | 3 -> " RR="
                  | _ -> assert false);
                Buffer.add_char buf (if b then 'Y' else 'N'))
              Parameters.allowed_pairing_orientations;
            Buffer.add_string buf " }";
            Buffer.contents buf),
        fun what ->
          let b =
            match what with
            | "--set-orientation" -> true
            | "--unset-orientation" -> false
            | _ -> assert false
          in
          match TA.get_parameter () with
          | "FF" -> Parameters.allowed_pairing_orientations.(0) <- b
          | "FR" -> Parameters.allowed_pairing_orientations.(1) <- b
          | "RF" -> Parameters.allowed_pairing_orientations.(2) <- b
          | "RR" -> Parameters.allowed_pairing_orientations.(3) <- b
          | w ->
              TA.usage ();
              Printf.sprintf
                "Expected pairing orientation FF|FR|RF|RR, found '%s'" w
              |> failwith );
      ( [ "-n"; "--maximum-pairing-number" ],
        Some "<positive_integer>",
        [
          "maximum number of pairings to be made.";
          "Ends will be output separately if threshold exceeded";
        ],
        TA.Default (fun () -> string_of_int Defaults.max_pairing_number),
        fun _ -> Parameters.max_pairing_number := TA.get_parameter_int_pos () );
      ( [ "--accept-unpairable-names" ],
        None,
        [
          "pair input even though the names of the ends do not appear";
          "to be consistent with usual conventions.";
          "The name of the first end will be used for the pair";
        ],
        TA.Default (fun () -> string_of_bool Defaults.accept_unpairable_names),
        fun _ -> Parameters.accept_unpairable_names := true );
      ( [ "-p"; "--unpaired-print-mode" ],
        Some "0|1|2",
        [
          "how to print unpaired reads.";
          "0: discard";
          "1: as one joint paired-end record without alignments";
          "2: as two separate single-end records";
        ],
        TA.Default
          (fun () ->
            AlignmentPairer.UnpairedPrintMode.to_string
              Defaults.unpaired_print_mode),
        fun _ ->
          Parameters.unpaired_print_mode :=
            AlignmentPairer.UnpairedPrintMode.of_string (TA.get_parameter ()) );
      ([], None, [ "Input/Output" ], TA.Optional, fun _ -> ());
      ( [ "--lines-per-block" ],
        Some "<positive_integer>",
        [ "number of lines to be processed per block" ],
        TA.Default (fun () -> string_of_int Defaults.lines_per_block),
        fun _ -> Parameters.lines_per_block := TA.get_parameter_int_pos () );
      ( [ "-i"; "--input" ],
        Some "<input_file>",
        [
          "name of GEM MAP input file.";
          "More than one can be specified;";
          "exactly two, or one interleaved, must be specified in pair mode";
        ],
        TA.Default (fun () -> "stdin"),
        fun _ -> TA.get_parameter () |> Tools.Misc.accum Parameters.inputs );
      ( [ "-o"; "--output" ],
        Some "<output_file>",
        [ "name of generated GEM MAP output file" ],
        TA.Default (fun () -> "stdout"),
        fun _ -> Parameters.output := TA.get_parameter () );
      ([], None, [ "Miscellaneous" ], TA.Optional, fun _ -> ());
      ( [ "--chunk-size" ],
        Some "<positive_integer>",
        [ "number of megabytes per chunk (64 can be a good value)";
        "If not providing `-i` only reads a single chunk"; ],
        TA.Mandatory,(* (fun () -> string_of_int Defaults.chunk_size), *)
        fun _ -> Parameters.chunk_size := TA.get_parameter_int_pos () );
      ( [ "-t"; "--threads" ],
        Some "<positive_integer>",
        [ "number of computing threads to be used" ],
        TA.Default (fun () -> string_of_int Defaults.threads),
        fun _ -> Parameters.threads := TA.get_parameter_int_pos () );
      ( [ "-v"; "--verbose" ],
        None,
        [ "set verbose execution" ],
        TA.Default (fun () -> string_of_bool Defaults.verbose),
        fun _ -> Parameters.verbose := true );
      ( [ "-h"; "--help" ],
        None,
        [ "print syntax and exit" ],
        TA.Optional,
        fun _ ->
          TA.usage ();
          exit 1 );
    ];
  let working_mode = Tools.Misc.unbox_opt !Parameters.working_mode
  and interleaved = ref false in
  (match working_mode with
  | Parameters.WorkingMode.Pair -> (
      match List.length !Parameters.inputs with
      | 0 | 1 -> interleaved := true (* 0 means stdin *)
      | 2 -> interleaved := false
      | w ->
          TA.usage ();
          Printf.sprintf
            "Expected 1 (interleaved) or 2 (matched) inputs, found %d" w
          |> failwith)
  | Merge -> ());
  let interleaved = !interleaved in
  let inputs =
    if !Parameters.inputs = [] then [| stdin |]
    else List.rev !Parameters.inputs |> List.map open_in |> Array.of_list
  and output =
    if !Parameters.output = "" then stdout else open_out !Parameters.output
  in
  let print_num_lines what =
    Printf.sprintf "%d %s %s" !Misc.input_line_num
      (Tools.Misc.pluralize_int "line" !Misc.input_line_num)
      what
    |> Tools.Misc.pteprintf "%s\n%!"
  in
  let lines_per_block = !Parameters.lines_per_block
  and eof = ref false
  and processing_buffer = Buffer.create 16777216 in
  let current_offset = ref 0 in
  Tools.Parallel.process_stream_chunkwise_with_lithops
    !Parameters.chunk_size
    (match !Parameters.inputs with [] -> None | first :: _ -> Some first)
    Tools.ElementDetector.map_se_re
    (fun input_string ->
      let read_until_next_newline () =
        try
          let index = String.index_from input_string !current_offset '\n' in
          let line_length = index - !current_offset in
          let read_line = String.sub input_string !current_offset line_length in
          current_offset := !current_offset + line_length + 1;
          incr Misc.input_line_num;
          read_line
        with _ -> raise End_of_file
      in

      if not !eof then (
        if !Misc.input_line_num mod lines_per_block = 0 then
          print_num_lines "read";
        let buf = ref [] and base_input_line_num = !Misc.input_line_num + 1 in
        (try
           if interleaved then
             for i = 1 to lines_per_block do
               let line_1 = read_until_next_newline () in
               let line_2 =
                 try read_until_next_newline ()
                 with End_of_file ->
                   (* When reading the second line we must always succeed *)
                   Misc.line__failwith
                     "Your interleaved file contains an odd number of lines, \
                      pairing aborted"
               in
               Tools.Misc.accum buf [| line_1; line_2 |]
             done
           else
             for i = 1 to lines_per_block do
               [| read_until_next_newline () |] |> Tools.Misc.accum buf
             done
         with End_of_file ->
           print_num_lines "read";
           eof := true);
        (base_input_line_num, !buf))
      else raise End_of_file)
    (fun (base_input_line_num, reads_buffer) ->
      (* We want to keep counters realiable in case of error *)
      Misc.input_line_num := base_input_line_num;
      let lines_processed = ref 0 in
      Buffer.clear processing_buffer;
      List.iter
        (fun lines ->
          (match working_mode with
          | Parameters.WorkingMode.Merge ->
              AlignmentMerger.merge ~verbose:!Parameters.verbose
                Misc.line__failwith processing_buffer lines
          | Pair ->
              AlignmentPairer.pair ~verbose:!Parameters.verbose
                ~accept_unpairable_names:!Parameters.accept_unpairable_names
                ~unpaired_print_mode:!Parameters.unpaired_print_mode
                !Parameters.max_pairing_distance
                Parameters.allowed_pairing_orientations
                !Parameters.max_pairing_number
                Misc.line__failwith processing_buffer lines);

          incr lines_processed)
        (List.rev reads_buffer);
      (!lines_processed, Buffer.contents processing_buffer))
    (fun (processed, buf) ->
      if !Misc.input_line_num = 0 then print_num_lines "reduced";
      let old_periods = !Misc.input_line_num / lines_per_block in
      Misc.input_line_num := !Misc.input_line_num + processed;
      if !Misc.input_line_num / lines_per_block > old_periods then
        print_num_lines "reduced";
      Printf.fprintf output "%s%!" buf)
    !Parameters.threads;
  ignore (print_num_lines "processed");
  (* Cleanup actions *)
  close_out output;
  Array.iter close_in inputs
