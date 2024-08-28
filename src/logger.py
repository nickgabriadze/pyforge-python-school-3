import logging
from os import getenv


def setupLogging():
    SERVER_ID = getenv("SERVER_ID", "1")

    logging.basicConfig(level=logging.INFO,
                        format=f'[DATE] %(asctime)s - [SERVER ID] {SERVER_ID} - [%(levelname)s] - %(message)s',
                        datefmt="%Y-%m-%d %H:%M:%S",
                        filename='fastapi.log')

    return logging
