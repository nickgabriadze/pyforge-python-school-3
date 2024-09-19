import os
from pydantic_settings import BaseSettings, SettingsConfigDict
from alembic.config import Config


class Settings(BaseSettings):
    DB_HOST: str
    DB_PORT: int
    DB_NAME: str
    DB_USER: str
    DB_PASSWORD: str
    model_config = SettingsConfigDict(
        extra='allow',
        env_file=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".env")
    )


settings = Settings()


def get_db_url():
    return (
        f"postgresql+asyncpg://{settings.DB_USER}:{settings.DB_PASSWORD}@"
        f"{settings.DB_HOST}:{settings.DB_PORT}/{settings.DB_NAME}"
    )


def get_alembic_config():
    config = Config("alembic.ini")
    config.set_main_option("sqlalchemy.url", settings.db_url)
    return config
