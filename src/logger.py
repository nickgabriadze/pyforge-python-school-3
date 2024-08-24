import logging


def setupLogging():
    logging.basicConfig(level=logging.INFO,
                        format='[DATE] %(asctime)s - [%(levelname)s] - %(message)s',
                        datefmt="%Y-%m-%d %H:%M:%S",
                        filename='fastapi.log')

    return logging
