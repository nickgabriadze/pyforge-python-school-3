from sqlalchemy.orm import DeclarativeBase, declared_attr, mapped_column
from sqlalchemy.ext.asyncio import create_async_engine, async_sessionmaker, AsyncAttrs
from typing import Annotated
from src.config import get_db_url

SQLALCHEMY_DB_URL = get_db_url()

engine = create_async_engine(SQLALCHEMY_DB_URL)

async_session_maker = async_sessionmaker(engine, expire_on_commit=False)

pubchem_pk = Annotated[str, mapped_column(primary_key=True)]
str_uniq = Annotated[str, mapped_column(unique=True, nullable=False)]
str_null_true = Annotated[str, mapped_column(nullable=True)]


class Base(AsyncAttrs, DeclarativeBase):
    __abstract__ = True

    @declared_attr.directive
    def __tablename__(cls) -> str:
        return f"{cls.__name__.lower()}s"
