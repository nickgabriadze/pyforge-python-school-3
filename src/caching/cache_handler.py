import json
import redis


def get_cached_result(client, key: str):
    result = client.get(key)

    if result:
        return json.loads(result)

    return None


def remove_cache(client, key:str):
    cacheExists = get_cached_result(client, key)
    if cacheExists:
        client.delete(key)

def set_cache(client, key: str, value: dict, exp: int = 60):
    client.setex(key, exp, json.dumps(value))

