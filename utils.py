import json
import functools
import configparser
from timeit import default_timer as timer
from colorama import init, Fore


def get_host_data() -> dict:
    """
        Reads Edwards' dataset host information.
    """
    with open("host.json", "r") as fh:
        host_data = json.load(fh)
    return host_data


def get_virus_data() -> dict:
    with open("virus.json", "r") as fh:
        virus_data = json.load(fh)
    return virus_data


def get_config() -> dict:
    """
        Reads important settings and variables to run the analysis properly
    """
    config = configparser.ConfigParser()
    config.read("config.cfg")
    config_dict = {section: dict(config.items(section)) for section in config.sections()}
    return config_dict


def get_config_obj():
    config = configparser.ConfigParser()
    config.read("config.cfg")
    return config


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

