import time

from src.celery_worker import celery

@celery.task(name="add_task")
def add_task(a: int, x: int, y: int):
    time.sleep(a)
    return x + y