import time
import asyncio
from src.celery_worker import celery


@celery.task(name="add_task")
def add_task(a: int, x: int, y: int):
    time.sleep(a)
    return x + y


@celery.task(name="sub_search_task")
def sub_search_task(mol: str, limit: int):
    from src.main import get_sub_match
    # Call the async function and get results
    loop = asyncio.get_event_loop()
    sub_matches = loop.run_until_complete(get_sub_match(mol, limit))
    # Convert the results to a dictionary or a proper JSON structure
    return sub_matches
