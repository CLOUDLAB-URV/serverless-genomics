(* This library is
    (C) 2015-2017 Paolo Ribeca, <paolo.ribeca@gmail.com>.
   It contains several general-purpose utilities, in particular:
    * a module to parse command-line options
    * a module to arbitrarily parallelize streams following a
       reader-workers-writer model.
   Please do not redistribute, or contact the author before doing so *)

module Misc : sig
  (* General utilities *)
  val round : float -> float

  (* *)
  val unbox_opt : 'a option -> 'a
  val unbox_opt_def : 'a -> 'a option -> 'a

  (* *)
  val accum : 'a list ref -> 'a -> unit
  val pop : 'a list ref -> 'a

  (* *)
  val pluralize : one:'a -> string -> 'a -> string
  val pluralize_int : string -> int -> string
  val pluralize_float : string -> float -> string

  (* *)
  type mode_t = Time | Space | Empty

  val tfprintf :
    ?mode:mode_t -> out_channel -> ('a, out_channel, unit) format -> 'a

  val tprintf : ?mode:mode_t -> ('a, out_channel, unit) format -> 'a
  val teprintf : ?mode:mode_t -> ('a, out_channel, unit) format -> 'a
  val pfprintf : out_channel -> ('a, out_channel, unit) format -> 'a
  val pprintf : ('a, out_channel, unit) format -> 'a
  val peprintf : ('a, out_channel, unit) format -> 'a

  val ptfprintf :
    ?mode:mode_t -> out_channel -> ('a, out_channel, unit) format -> 'a

  val ptprintf : ?mode:mode_t -> ('a, out_channel, unit) format -> 'a
  val pteprintf : ?mode:mode_t -> ('a, out_channel, unit) format -> 'a
end = struct
  (* General utilities *)
  let round x = if x >= 0. then floor (x +. 0.5) else ceil (x -. 0.5)

  (* *)
  let unbox_opt = function None -> assert false | Some r -> r
  let unbox_opt_def def = function None -> def | Some r -> r

  (* *)
  let accum rl el = rl := el :: !rl

  let pop rl =
    match !rl with
    | [] -> raise Not_found
    | hd :: tail ->
        rl := tail;
        hd

  (* *)
  let pluralize (type a) ~(one : a) s n = if n = one then s else s ^ "s"
  let pluralize_int = pluralize ~one:1
  let pluralize_float = pluralize ~one:1.

  (* *)
  type mode_t = Time | Space | Empty

  let tfprintf ?(mode = Time) ch =
    let t = Unix.localtime (Unix.time ()) in
    (match mode with
    | Time ->
        Printf.fprintf ch "%s %s %2d %02d:%02d:%02d %4d -- "
          (match t.Unix.tm_wday with
          | 0 -> "Sun"
          | 1 -> "Mon"
          | 2 -> "Tue"
          | 3 -> "Wed"
          | 4 -> "Thu"
          | 5 -> "Fri"
          | 6 -> "Sat"
          | _ -> assert false)
          (match t.Unix.tm_mon with
          | 0 -> "Jan"
          | 1 -> "Feb"
          | 2 -> "Mar"
          | 3 -> "Apr"
          | 4 -> "May"
          | 5 -> "Jun"
          | 6 -> "Jul"
          | 7 -> "Aug"
          | 8 -> "Sep"
          | 9 -> "Oct"
          | 10 -> "Nov"
          | 11 -> "Dec"
          | _ -> assert false)
          t.Unix.tm_mday t.Unix.tm_hour t.Unix.tm_min t.Unix.tm_sec
          (1900 + t.Unix.tm_year)
    | Space -> Printf.fprintf ch "                         -- "
    | Empty -> ());
    Printf.fprintf ch

  let tprintf ?(mode = Time) = tfprintf ~mode stdout
  let teprintf ?(mode = Time) = tfprintf ~mode stderr

  let proto_pfprintf ?(mode = Time) f ch =
    (match mode with
    | Time -> Unix.getpid () |> Printf.fprintf ch "[%07d] "
    | Space -> Printf.fprintf ch "          "
    | Empty -> ());
    f ch

  let pfprintf ch = proto_pfprintf ~mode:Time Printf.fprintf ch
  let pprintf w = pfprintf stdout w
  let peprintf w = pfprintf stderr w
  let ptfprintf ?(mode = Time) = proto_pfprintf ~mode (tfprintf ~mode)
  let ptprintf ?(mode = Time) = ptfprintf ~mode stdout
  let pteprintf ?(mode = Time) = ptfprintf ~mode stderr
end

module Split = struct
  let tab_re = Str.regexp "\t"
  and space_re = Str.regexp " "
  and comma_re = Str.regexp ","
  and colon_re = Str.regexp ":"
  and newline_re = Str.regexp "\n"
  and underscore_re = Str.regexp "_"
  and fasta_name_re = Str.regexp "^>"

  let as_list = Str.split
  let as_array re s = Array.of_list (Str.split re s)
end

module Argv : sig
  type class_t =
    | Mandatory (* Implies no default *)
    | Optional (* Implies no default *)
    | Default of (unit -> string)
  (* Implies optional - the function prints the default *)

  (* The specs from which usage is produced are a tuple with the following elements:
      (1) equivalent option names
      (2) optional explanation for the argument(s)
      (3) help text lines
      (4) mandatory/optional * default/no-default class
      (5) parsing action when option is encountered *)
  type spec_t =
    string list * string option * string list * class_t * (string -> unit)

  val usage : ?output:out_channel -> unit -> unit
  val get_parameter : unit -> string
  val get_parameter_boolean : unit -> bool
  val get_parameter_int : unit -> int
  val get_parameter_float : unit -> float
  val get_parameter_int_pos : unit -> int
  val get_parameter_float_pos : unit -> float
  val get_parameter_int_non_neg : unit -> int
  val get_parameter_float_non_neg : unit -> float
  val get_parameter_float_fraction : unit -> float

  (* Consumes and returns all the parameters which are left on the command line *)
  val get_remaining_parameters : unit -> string array
  val parse : spec_t list -> unit
end = struct
  type class_t = Mandatory | Optional | Default of (unit -> string)

  type spec_t =
    string list * string option * string list * class_t * (string -> unit)

  let argv = Sys.argv
  let i = ref 1
  let buf = ref ("Usage:\n  " ^ argv.(0) ^ "\n")
  let usage ?(output = stdout) () = Printf.fprintf output "%s" !buf

  let error ?(output = stderr) f_n msg =
    usage ~output ();
    Printf.fprintf output "Tools.Argv.%s: %s\n%!" f_n msg;
    exit 1

  let template_get n what f () =
    try f ()
    with _ ->
      error n ("Option '" ^ argv.(!i - 1) ^ "' needs a " ^ what ^ " parameter")

  let template_filter f g () =
    let res = f () in
    if res |> g then res else raise Not_found

  let get_parameter =
    template_get "get_parameter" "" (fun () ->
        incr i;
        argv.(!i))

  let get_parameter_boolean =
    template_get "get_parameter_boolean" "boolean" (fun () ->
        get_parameter () |> bool_of_string)

  let get_parameter_int =
    template_get "get_parameter_int" "integer" (fun () ->
        get_parameter () |> int_of_string)

  let get_parameter_float =
    template_get "get_parameter_float" "float" (fun () ->
        get_parameter () |> float_of_string)

  let get_parameter_int_pos =
    template_get "get_parameter_int_pos" "positive integer"
      (template_filter get_parameter_int (fun x -> x > 0))

  let get_parameter_float_pos =
    template_get "get_parameter_float_pos" "positive float"
      (template_filter get_parameter_float (fun x -> x > 0.))

  let get_parameter_int_non_neg =
    template_get "get_parameter_int_non_neg" "non-negative integer"
      (template_filter get_parameter_int (fun x -> x >= 0))

  let get_parameter_float_non_neg =
    template_get "get_parameter_float_non_neg" "non-negative float"
      (template_filter get_parameter_float (fun x -> x >= 0.))

  let get_parameter_float_fraction =
    template_get "get_parameter_float_fraction" "float between 0 and 1"
      (template_filter get_parameter_float (fun x -> x >= 0. && x <= 1.))

  let get_remaining_parameters () =
    let len = Array.length argv in
    let res = Array.sub argv (!i + 1) (len - !i - 1) in
    i := len;
    res

  let parse specs =
    let module StringMap = Map.Make (String) in
    let module StringSet = Set.Make (String) in
    let table = ref StringMap.empty and mandatory = ref StringSet.empty in
    List.iteri
      (fun i (opts, vl, help, clss, act) ->
        buf := !buf ^ if opts <> [] then "  " else " ";
        if opts = [] && help = [] then
          error "parse" ("Malformed initializer for option #" ^ string_of_int i);
        List.iteri
          (fun i opt ->
            if i > 0 then buf := !buf ^ "|";
            buf := !buf ^ opt;
            if StringMap.mem opt !table then
              error "parse"
                ("Duplicate command line option '" ^ opt ^ "' in table");
            if clss = Mandatory then (
              let repr =
                List.fold_left
                  (fun a b -> a ^ (if a = "" then "" else "|") ^ b)
                  "" opts
              in
              mandatory := StringSet.add repr !mandatory;
              table :=
                StringMap.add opt
                  (fun arg ->
                    mandatory := StringSet.remove repr !mandatory;
                    act arg)
                  !table)
            else table := StringMap.add opt act !table)
          opts;
        buf :=
          !buf
          ^ (match vl with None -> "" | Some vl -> " " ^ vl)
          ^ if opts <> [] then "\n" else "";
        List.iter
          (if opts <> [] then fun help -> buf := !buf ^ "     " ^ help ^ "\n"
          else fun help -> buf := !buf ^ help ^ "\n")
          help;
        match clss with
        | Mandatory -> buf := !buf ^ "     (mandatory)\n"
        | Optional -> ()
        | Default def -> buf := !buf ^ "     (default='" ^ def () ^ "')\n")
      specs;
    let len = Array.length argv in
    while !i < len do
      let arg = argv.(!i) in
      (try StringMap.find arg !table
       with Not_found -> error "parse" ("Unknown option '" ^ argv.(!i) ^ "'"))
        arg;
      incr i
    done;
    if !mandatory <> StringSet.empty then
      StringSet.iter
        (fun opt -> error "parse" ("Option '" ^ opt ^ "' is mandatory"))
        !mandatory
end

module ElementDetector = struct
  (* Better to try to catch the overall structure of the element *)
  (** WARNING: Matching a single char is not programmed properly *)
  let fastq_re = Str.regexp "^@.*\n*.+\n\\+.*\n.+$"

  (* TODO: improve *)
  and map_se_re = Str.regexp "^[A-Z0-9]+\\.[0-9].*$"
end

module Parallel : sig
  val process_stream_chunkwise :
    ?buffered_chunks_per_thread:int ->
    (* Beware: for everything to terminate properly, f shall raise End_of_file when done.
       Side effects are propagated within f (not exported) and within h (exported) *)
    (unit -> 'a) ->
    ('a -> 'b) ->
    ('b -> unit) ->
    int ->
    unit

  val process_stream_chunkwise_with_lithops :
    int ->
    string option ->
    Str.regexp ->
    (string -> 'a) ->
    ('a -> 'b) ->
    ('b -> unit) ->
    int ->
    unit

  val process_stream_linewise :
    ?buffered_chunks_per_thread:int ->
    ?max_memory:int ->
    ?string_buffer_memory:int ->
    ?input_line:(in_channel -> string) ->
    ?verbose:bool ->
    in_channel ->
    (Buffer.t -> int -> string -> unit) ->
    out_channel ->
    int ->
    unit
end = struct
  module Int_map = Map.Make (struct
    type t = int

    let compare = compare
  end)

  let process_stream_chunkwise ?(buffered_chunks_per_thread = 10)
      (f : unit -> 'a) (g : 'a -> 'b) (h : 'b -> unit) threads =
    if buffered_chunks_per_thread < 1 || threads < 1 then
      raise (Failure "parallel_process_stream_chunkwise");

    let red_threads = threads - 1 in
    (* I am the ouptut process *)
    let close_pipe (pipe_in, pipe_out) =
      Unix.close pipe_in;
      Unix.close pipe_out
    in
    let close_pipes_in =
      Array.iter (fun (pipe_in, pipe_out) -> Unix.close pipe_in)
    and close_pipes_out =
      Array.iter (fun (pipe_in, pipe_out) -> Unix.close pipe_out)
    and close_pipes = Array.iter close_pipe
    and get_stuff_for_select pipes =
      let pipes = Array.map fst pipes in
      ( Array.to_list pipes,
        let dict = Hashtbl.create (Array.length pipes) in
        Array.iteri
          (fun i pipe ->
            assert (not (Hashtbl.mem dict pipe));
            Hashtbl.add dict pipe i)
          pipes;
        dict )
    and w_2_o_pipes = Array.init threads (fun _ -> Unix.pipe ())
    and o_2_w_pipes = Array.init threads (fun _ -> Unix.pipe ()) in
    flush_all ();
    match Unix.fork () with
    | 0 ->
        (* Child *)
        (* I am the input process *)
        let i_2_w_pipes = Array.init threads (fun _ -> Unix.pipe ())
        and w_2_i_pipes = Array.init threads (fun _ -> Unix.pipe ()) in
        for i = 0 to red_threads do
          match Unix.fork () with
          | 0 ->
              (* Child *)
              (* I am a worker.
                 I only keep my own pipes open *)
              let i_2_w_pipe_in, w_2_i_pipe_out, o_2_w_pipe_in, w_2_o_pipe_out =
                let i_2_w_pipe_in = ref Unix.stdin
                and w_2_i_pipe_out = ref Unix.stdout
                and o_2_w_pipe_in = ref Unix.stdin
                and w_2_o_pipe_out = ref Unix.stdout in
                for ii = 0 to red_threads do
                  if ii = i then (
                    let pipe_in, pipe_out = i_2_w_pipes.(ii) in
                    i_2_w_pipe_in := pipe_in;
                    Unix.close pipe_out;
                    let pipe_in, pipe_out = w_2_i_pipes.(ii) in
                    Unix.close pipe_in;
                    w_2_i_pipe_out := pipe_out;
                    let pipe_in, pipe_out = o_2_w_pipes.(ii) in
                    o_2_w_pipe_in := pipe_in;
                    Unix.close pipe_out;
                    let pipe_in, pipe_out = w_2_o_pipes.(ii) in
                    Unix.close pipe_in;
                    w_2_o_pipe_out := pipe_out)
                  else (
                    close_pipe i_2_w_pipes.(ii);
                    close_pipe w_2_i_pipes.(ii);
                    close_pipe o_2_w_pipes.(ii);
                    close_pipe w_2_o_pipes.(ii))
                done;
                ( !i_2_w_pipe_in,
                  !w_2_i_pipe_out,
                  !o_2_w_pipe_in,
                  !w_2_o_pipe_out )
              in
              let i_2_w = Unix.in_channel_of_descr i_2_w_pipe_in
              and w_2_i = Unix.out_channel_of_descr w_2_i_pipe_out
              and o_2_w = Unix.in_channel_of_descr o_2_w_pipe_in
              and w_2_o = Unix.out_channel_of_descr w_2_o_pipe_out in
              (* My protocol is:
                 (1) process a chunk more from the input
                 (2) notify the output process that a result is ready
                 (3) when the output process asks for it, post the result *)
              let probe_output () =
                (*ignore (Unix.select [o_2_w_pipe_in] [] [] (-1.));*)
                ignore (input_byte o_2_w)
              and initial = ref true in
              while true do
                (* Try to get one more chunk.
                   Signal the input process that I am idle *)
                output_byte w_2_i 0;
                flush w_2_i;
                (* Get & process a chunk *)
                match input_byte i_2_w with
                | 0 ->
                    (* EOF reached *)
                    (* Did the output process ask for a notification? *)
                    if not !initial then
                      (* The first time, we notify anyway to avoid crashes *)
                      probe_output ();
                    (* Notify that EOF has been reached *)
                    output_binary_int w_2_o (-1);
                    flush w_2_o;
                    (* Did the output process switch me off? *)
                    probe_output ();
                    (* Commit suicide *)
                    Unix.close i_2_w_pipe_in;
                    Unix.close w_2_i_pipe_out;
                    Unix.close o_2_w_pipe_in;
                    Unix.close w_2_o_pipe_out;
                    exit 0
                | 1 ->
                    (* OK, one more token available *)
                    (* Get the chunk *)
                    let (chunk_id, data) = (input_value i_2_w : int * 'a) in
                    (* Process the chunk *)
                    let data = g data in
                    (* Did the output process ask for a notification? *)
                    if not !initial then
                      (* The first time, we notify anyway to avoid crashes *)
                      probe_output ()
                    else initial := false;
                    (* Tell the output process what we have *)
                    output_binary_int w_2_o chunk_id;
                    flush w_2_o;
                    (* Did the output process request data? *)
                    probe_output ();
                    (* Send the data to output *)
                    output_value w_2_o data;
                    flush w_2_o
                | _ -> assert false
              done
          | _ -> () (* Parent *)
        done;
        (* I am the input process.
           I do not care about output process pipes *)
        close_pipes w_2_o_pipes;
        close_pipes o_2_w_pipes;
        close_pipes_in i_2_w_pipes;
        close_pipes_out w_2_i_pipes;
        let w_2_i_pipes_for_select, w_2_i_dict =
          get_stuff_for_select w_2_i_pipes
        and w_2_i =
          Array.map
            (fun (pipe_in, _) -> Unix.in_channel_of_descr pipe_in)
            w_2_i_pipes
        and i_2_w =
          Array.map
            (fun (_, pipe_out) -> Unix.out_channel_of_descr pipe_out)
            i_2_w_pipes
        in
        (* My protocol is:
           (1) read a chunk
           (2) read a thread id
           (3) post the chunk to the correspondng pipe. The worker will consume it *)
        let chunk_id = ref 0 and off = ref 0 in
        while !off < threads do
          let ready, _, _ = Unix.select w_2_i_pipes_for_select [] [] (-1.) in
          List.iter
            (fun ready ->
              let w_id = Hashtbl.find w_2_i_dict ready in
              ignore (input_byte w_2_i.(w_id));
              let i_2_w = i_2_w.(w_id) in
              try
                if !off > 0 then raise End_of_file;
                let payload = f () in
                output_byte i_2_w 1;
                (* OK to transmit, we have not reached EOF yet *)
                flush i_2_w;
                output_value i_2_w (!chunk_id, payload);
                flush i_2_w;
                incr chunk_id
              with End_of_file ->
                output_byte i_2_w 0;
                (* Nothing to transmit *)
                flush i_2_w;
                incr off)
            ready
        done;
        (* Waiting to be switched off *)
        ignore (Unix.select [ fst w_2_i_pipes.(0) ] [] [] (-1.));
        close_pipes_out i_2_w_pipes;
        close_pipes_in w_2_i_pipes;
        exit 0
    | _ ->
        (* I am the output process *)
        close_pipes_in o_2_w_pipes;
        close_pipes_out w_2_o_pipes;
        let w_2_o_pipes_for_select, w_2_o_dict =
          get_stuff_for_select w_2_o_pipes
        and w_2_o =
          Array.map
            (fun (pipe_in, _) -> Unix.in_channel_of_descr pipe_in)
            w_2_o_pipes
        and o_2_w =
          Array.map
            (fun (_, pipe_out) -> Unix.out_channel_of_descr pipe_out)
            o_2_w_pipes
        and buffered_chunks = buffered_chunks_per_thread * threads
        and next = ref 0
        and queue = ref Int_map.empty
        and buf = ref Int_map.empty
        and off = ref 0 in
        let o_2_stdout_pipe_out = ref Unix.stdout in
        let o_2_stdout = Unix.out_channel_of_descr !o_2_stdout_pipe_out in
        while !off < threads do
          (* Harvest new notifications *)
          let ready, _, _ = Unix.select w_2_o_pipes_for_select [] [] (-1.) in
          List.iter
            (fun ready ->
              let w_id = Hashtbl.find w_2_o_dict ready in
              let chunk_id = input_binary_int w_2_o.(w_id) in
              if chunk_id = -1 then (* EOF has been reached *)
                incr off
              else if not (Int_map.mem chunk_id !queue) then
                queue := Int_map.add chunk_id w_id !queue
              else assert (w_id = Int_map.find chunk_id !queue))
            ready;
          (* Fill the buffer *)
          let available = ref (buffered_chunks - Int_map.cardinal !buf) in
          assert (!available >= 0);
          (* If the needed chunk is there, we always fetch it *)
          if !queue <> Int_map.empty && fst (Int_map.min_binding !queue) = !next
          then incr available;
          while !available > 0 && !queue <> Int_map.empty do
            let chunk_id, w_id = Int_map.min_binding !queue in
            (* Tell the worker to send data *)
            output_byte o_2_w.(w_id) 0;
            flush o_2_w.(w_id);
            assert (not (Int_map.mem chunk_id !buf));
            buf := Int_map.add chunk_id (input_value w_2_o.(w_id) : 'b) !buf;
            (* Tell the worker to send the next notification *)
            output_byte o_2_w.(w_id) 0;
            flush o_2_w.(w_id);
            queue := Int_map.remove chunk_id !queue;
            decr available
          done;
          (* Output at most as many chunks at the number of workers *)
          available := threads;
          while !available > 0 && !buf <> Int_map.empty do
            let chunk_id, data = Int_map.min_binding !buf in
            if chunk_id = !next then (
              (* Worker role: output the received chunks from the workers directly to stdout *)
              if Sys.getenv_opt "DIFF" <> None then (
                output_value o_2_stdout data;
                flush o_2_stdout)
              else (* Worker role not set: reduce directly *)
                h data;

              buf := Int_map.remove chunk_id !buf;
              incr next;
              decr available)
            else available := 0 (* Force exit from the cycle *)
          done
        done;
        (* There might be chunks left in the buffer *)
        while !buf <> Int_map.empty do
          let chunk_id, data = Int_map.min_binding !buf in
          assert (chunk_id = !next);

          (* Worker role: output the received chunks from the workers directly to stdout *)
          if Sys.getenv_opt "DIFF" <> None then (
            output_value o_2_stdout data;
            flush o_2_stdout)
          else (* Worker role not set: reduce directly *)
            h data;

          buf := Int_map.remove chunk_id !buf;
          incr next
        done;
        (* Switch off all the workers *)
        let red_threads = threads - 1 in
        for ii = 0 to red_threads do
          output_byte o_2_w.(ii) 0;
          flush o_2_w.(ii)
        done;
        close_pipes_out o_2_w_pipes;
        close_pipes_in w_2_o_pipes

  let partition_input_file filename chunk_threshold obj_chunk_size =
    let stat = Unix.stat filename in
    let start_offset = ref 0 and ranges = ref [] in
    (* create ranges with fixed size except for the last one *)
    while !start_offset < stat.st_size - obj_chunk_size do
      let range =
        (!start_offset, !start_offset + obj_chunk_size + chunk_threshold)
      in
      ranges := range :: !ranges;
      start_offset := !start_offset + obj_chunk_size
    done;
    (* add the last exact range to cover until the end of the file *)
    ranges := (!start_offset, stat.st_size) :: !ranges;
    List.rev !ranges

  let byte_range_stream_from_file filename ranges =
    let file_chan = open_in_bin filename in
    Stream.from (fun range_index ->
        if List.length ranges <= range_index then None
        else
          let start_offset, end_offset = List.nth ranges range_index in
          seek_in file_chan start_offset;
          Misc.teprintf "Range start %d, end %d\n" start_offset end_offset;
          let chunk_size = end_offset - start_offset in
          Some (really_input_string file_chan chunk_size))

  let heal_chunk_with_regexp
      ((buf : string), (chunk_size : int), (element_delimiter : Str.regexp)) =
    (* Lithops adds the last byte from the previous chunk in the first byte in the buffer *)
    (* except for the first chunk (which has no previous chunk) *)
    (* However, right now, matching single chars is NOT supported. *)
    let first_elem_start =
      try Str.search_forward element_delimiter buf 0
      with Not_found ->
        Misc.teprintf "Errored string start:%s...\n%!" (Str.first_chars buf 250);
        invalid_arg "Impossible to heal the function. First element not found."
    in
    let string_buf_size = String.length buf in
    if string_buf_size >= chunk_size then
      let last_elem_end =
        try Str.search_forward element_delimiter buf chunk_size
        with Not_found ->
          Misc.teprintf "Errored string end + 100 threshold chars:...%s\n%!"
            (String.sub buf first_elem_start
               (string_buf_size - first_elem_start + 500));
          invalid_arg "Impossible to heal the function. Last element not found."
      in
      String.sub buf first_elem_start (last_elem_end - first_elem_start)
    else
      (* If the last chunk is smaller than the others just return until the end *)
      String.sub buf first_elem_start (string_buf_size - first_elem_start)

  let process_stream_chunkwise_with_lithops chunk_size_in_mb
      (filename : string option) (elem_detector : Str.regexp)
      (chunker : string -> 'a) (worker : 'a -> 'b) (reducer : 'b -> unit)
      threads =
    let chunk_size = chunk_size_in_mb * 1024 * 1024 in
    (* The threshold is a lithops constant *)
    let chunk_threshold = 128 * 1024 in

    if Option.is_some filename then
      (* Lithops is not running. We want to partition a *file* simulating lithops. *)
      let filename = Option.get filename in
      let ranges = partition_input_file filename chunk_threshold chunk_size in
      let stream = byte_range_stream_from_file filename ranges in
      let healed string_buf =
        heal_chunk_with_regexp (string_buf, chunk_size, elem_detector)
      in
      let execute_par string_buf =
        process_stream_chunkwise
          (fun () -> chunker (healed string_buf))
          worker reducer threads
      in
      Stream.iter execute_par stream
    else
      (* Filename not set, read from stdin only. *)
      match Sys.getenv_opt "DIFF" with
      | Some "reducer" -> (
          (* Reducer is single threaded *)
          let w_2_o_pipe_in = ref Unix.stdin in
          let w_2_o = Unix.in_channel_of_descr !w_2_o_pipe_in in
          try
            while true do
              let data = (input_value w_2_o : 'b) in
              reducer data
            done
          with End_of_file -> ())
      | None | Some _ ->
          (* Assume receiving a single chunk *)
          let read_stdin_into_buf buffer_size =
            let input_buffer = Buffer.create buffer_size in
            try
              Buffer.add_channel input_buffer stdin buffer_size;
              input_buffer
            with End_of_file -> input_buffer
          in
          let read_stdin = read_stdin_into_buf (chunk_size + chunk_threshold) in
          let healed =
            heal_chunk_with_regexp
              (Buffer.contents read_stdin, chunk_size, elem_detector)
          in
          process_stream_chunkwise
            (fun () -> chunker healed)
            worker reducer threads

  let process_stream_linewise ?(buffered_chunks_per_thread = 10)
      ?(max_memory = 1000000000) ?(string_buffer_memory = 16777216)
      ?(input_line = input_line) ?(verbose = true) input
      (f : Buffer.t -> int -> string -> unit) output threads =
    let max_block_bytes = max_memory / (buffered_chunks_per_thread * threads) in
    (* Parallel section *)
    let read = ref 0
    and eof_reached = ref false
    and processing_buffer = Buffer.create string_buffer_memory
    and processed = ref 0
    and written = ref 0 in
    if verbose then Misc.teprintf "0 lines read\n";
    process_stream_chunkwise ~buffered_chunks_per_thread
      (fun () ->
        if not !eof_reached then (
          let bytes = ref 0 and read_base = !read and buf = ref [] in
          (try
             while !bytes < max_block_bytes do
               (* We read at least one line *)
               let line = input_line input in
               Misc.accum buf line;
               bytes := String.length line + !bytes;
               incr read
             done
           with End_of_file -> eof_reached := true);
          if verbose then
            Misc.teprintf "%d %s read\n" !read (Misc.pluralize_int "line" !read);
          (read_base, !read - read_base, !buf))
        else raise End_of_file)
      (fun (lines_base, lines, buf) ->
        Buffer.clear processing_buffer;
        processed := 0;
        List.iter
          (fun line ->
            f processing_buffer (lines_base + !processed) line;
            incr processed)
          (List.rev buf);
        assert (!processed = lines);
        if verbose then
          Misc.teprintf "%d more %s processed\n" lines
            (Misc.pluralize_int "line" lines);
        (lines_base, lines, Buffer.contents processing_buffer))
      (fun (_, buf_len, buf) ->
        written := !written + buf_len;
        Printf.fprintf output "%s%!" buf;
        if verbose then
          Misc.teprintf "%d %s written\n" !written
            (Misc.pluralize_int "line" !written))
      threads;
    flush output;
    if verbose then
      Misc.teprintf "%d %s out\n" !written (Misc.pluralize_int "line" !written)
end
