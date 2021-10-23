import json
import functools
from timeit import default_timer as timer
from colorama import init, Fore


def get_host_data():
    """
        Reads Edwards' dataset host information.
    """
    with open("host.json", "r") as fh:
        host_data = json.load(fh)
    return host_data


def time_this(func):
    """
        Decorator which returns information about execution of decorated function.
    """
    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        start = timer()
        values = func(*args, **kwargs)
        end = timer()
        runtime = end - start
        if values is None:
            print(f"{Fore.RED}{func.__name__!r} execution error")
        else:
            print(f"{Fore.GREEN}{func.__name__!r} executed successfully in {runtime:.6f} seconds")
            return values[0]
        return wrapper_timer

