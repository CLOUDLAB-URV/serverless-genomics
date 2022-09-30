from enum import Enum

class Prices(float, Enum):
    t2medium_second = 0.00001288888
    gbxms_price = 0.0000000167
    s3_price_per_operation = 0.000005 

class PriceEstimator:
    def __init__(self,fexec):
        self.fexec = fexec
        self.arr = []

    #Calculates the cost of lambda functions based on gb/s

    def lambda_calc(self,rtm):
        stats = [f.stats for f in self.fexec.futures]
        sum_total_time = sum([stat['worker_exec_time'] for stat in stats]) * 1000
        total_lambda_price = Prices.gbxms_price.value * sum_total_time * rtm  # Price GB/ms * sum of times in ms * 1 GB
        return total_lambda_price
    
    #Calculates the cost of PUT, COPY, POST, or LIST requests

    def s3_calc(self):
        i = 0
        s3_price = 0
        while i < len(self.arr):
            if i != 0:
                print(self.arr[i][0] - self.arr[i-1][0])
                s3_price = s3_price + (self.arr[i][0] - self.arr[i-1][0]) * self.arr[i][1] * Prices.s3_price_per_operation.value
            else:
                print(self.arr[i][0])
                s3_price = s3_price + (self.arr[i][0]) * self.arr[i][1] * Prices.s3_price_per_operation.value
            i=i+1
        return s3_price

    def s3_calc_multiple_fexec(self, num_op):
        s3_price=0
        s3_price = s3_price + (len(self.fexec.futures)) * num_op * Prices.s3_price_per_operation.value
        return s3_price
    
    #Calculates ec2 costs based on uptime
    def ec2_calc(self,execution_time):
        return Prices.t2medium_second.value * execution_time

    