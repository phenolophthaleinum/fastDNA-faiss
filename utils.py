import json
import functools
import configparser
import re
from timeit import default_timer as timer
from colorama import init, Fore


def get_host_data() -> dict:
    """Reads Edwards' dataset of host information.

    Returns:
        dict: Returns dictionary of data parsed from `host.json`
    """
    with open("host.json", "r") as fh:
        host_data = json.load(fh)
    return host_data


def get_virus_data() -> dict:
    """Reads Edwards' dataset of virus information.

    Returns:
        dict: Returns dictionary of data parsed from `virus.json`
    """
    with open("virus.json", "r") as fh:
        virus_data = json.load(fh)
    return virus_data


def get_config() -> dict:
    """Reads important settings and variables to run the analysis properly

    Returns:
        dict: Returns dictionary of program variables parsed from `config.cfg`
    """
    config = configparser.ConfigParser()
    config.read("config.cfg")
    config_dict = {section: dict(config.items(section)) for section in config.sections()}
    return config_dict


def get_config_obj() -> configparser.ConfigParser:
    """Creates program configuration object from configuration file

    Returns:
        ConfigParser: Returns ConfigParser object created from `config.cfg`
    """
    config = configparser.ConfigParser()
    config.read("config.cfg")
    return config


def make_tax_json():
    """Creates `tax.json` file which is a modified Edwards' `host.json` file where keys are represented as taxid but
    ncbi_id information is still retained as one of the values

    """
    host_data = get_host_data()
    keys = list(host_data.keys())
    for key in keys:
        taxid = host_data[key]["taxid"]
        host_data[key]["ncbi_id"] = host_data[key].pop("taxid")
        host_data[key]["ncbi_id"] = key
        data = host_data[key]
        host_data[taxid] = host_data.pop(key)
    with open("tax.json", "w") as fh:
        json.dump(host_data, fh, indent=4)


def get_tax_data() -> dict:
    """Reads modified dataset of host information (taxid as a key).

    Returns:
        dict: Returns dictionary of data parsed from `tax.json`
    """
    with open("tax.json", "r") as fh:
        tax_data = json.load(fh)
    return tax_data


def make_hostvir_json():
    """Creates `hostvir.json` file which is a modified Edwards' `host.json` file with additional information about
    infecting viruses. Each record has a value with a list of infecting viruses.

    """
    host_data = get_host_data()
    virus_data = get_virus_data()
    for host in host_data:
        lineage_name = host_data[host]['lineage_names'][6]
        viruses = []
        for virus in virus_data:
            if virus_data[virus]['host']['organism_name'] == lineage_name:
                viruses.append(virus)
            host_data[host]['virus_id'] = viruses
    with open('hostvir.json', 'w') as fh:
        json.dump(host_data, fh)


def get_hostvir_data() -> dict:
    """Reads modified dataset of host information (additional information about infecting viruses).

    Returns:
        dict: Returns dictionary of data parsed from `hostvir.json`
    """
    with open('hostvir.json', 'r') as fh:
        hostvir_data = json.load(fh)
    return hostvir_data


def make_hostname_json():
    """Creates `hostname.json` file which is a modified Edwards' `host.json` file where keys are represented as organism name but
    ncbi_id information is still retained as one of the values

    """
    host_data = get_host_data()
    keys = list(host_data.keys())
    for key in keys:
        temp_name = host_data[key]['lineage_names'][-1]
        host_data[key]["ncbi_id"] = key
        name = "_".join(re.split(' |\; |\. |\, ', temp_name))
        host_data[name] = host_data.pop(key)
    with open("hostname.json", "w") as fh:
        json.dump(host_data, fh, indent=4)


def get_hostname_data() -> dict:
    """Reads modified dataset of host information (organism name as a key).

    Returns:
        dict: Returns dictionary of data parsed from `hostname.json`
    """
    with open("hostname.json", "r") as fh:
        hostname_data = json.load(fh)
    return hostname_data


def time_this(func):
    """*UNUSED* Decorator which returns information about execution of decorated function.
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

