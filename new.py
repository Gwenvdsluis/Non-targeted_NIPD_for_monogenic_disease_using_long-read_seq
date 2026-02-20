#!/usr/bin/env/ pyhton3
import string
n=1
result = ""
while n > 0:
    n, remainder = divmod(n - 1, 26)
    result = string.ascii_uppercase[remainder] + result
    print(result)



