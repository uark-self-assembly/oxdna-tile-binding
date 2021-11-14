from .env_vars import *

def set_debug_level(level):
    DEBUG_LEVEL = level
    print(f'DEBUG_LEVEL set to {DEBUG_LEVEL}')

def get_debug_level():
    return DEBUG_LEVEL

def print_message(level, msg):
    if DEBUG_LEVEL >= level:
        print(msg)
        with open("debug-output.log", "a") as f:
            f.write(msg + "\n")
