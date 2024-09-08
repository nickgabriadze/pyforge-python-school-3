from sqlalchemy import asc, cast, Integer, func, delete
from sqlalchemy.future import select
from ..molecules.models import Molecule
from .database import async_session_maker


class MoleculesDAO:
    model = Molecule

    @classmethod
    async def get_all_molecules(cls, limit=100):
        async with async_session_maker() as session:
            # since I have PUBCHEM{number} as a primary column in my database, I have to use a different
            # approach for the ascending style
            query = select(cls.model).order_by(
                asc(
                    cast(func.substring(cls.model.pubchem_id, r'[0-9]+'), Integer)
                )
            ).limit(limit)
            molecules = await session.execute(query)
            return molecules.scalars().all()

    @classmethod
    async def add_smiles(cls, **molecule):
        async with async_session_maker() as session:
            async with session.begin():
                check_if_exists = select(cls.model).filter_by(smiles=molecule.get('smiles'))
                request = await session.execute(check_if_exists)
                data = request.scalar_one_or_none()
                if data is None:
                    new_smiles = cls.model(**molecule)
                    session.add(new_smiles)
                    await session.flush()
                    await session.commit()
                    return molecule.get("pubchem_id")
                else:
                    return None

    @classmethod
    async def update_molecule(cls, pubchem_id: str, new_mol_smiles: str):
        async with async_session_maker() as session:
            async with session.begin():
                query = select(cls.model).filter_by(pubchem_id=pubchem_id)
                query_if_exists = select(cls.model).filter_by(smiles=new_mol_smiles)
                result = await session.execute(query)
                result_if_exists = await session.execute(query_if_exists)
                molecule_to_update = result.scalar_one_or_none()
                molecule_if_exists = result_if_exists.one_or_none()

                if molecule_to_update is None:
                    return None
                elif molecule_if_exists is not None:
                    return None

                molecule_to_update.smiles = new_mol_smiles
                await session.commit()
                return new_mol_smiles

    @classmethod
    async def get_molecule_by_pubchem_id(cls, pubchem_id: str):
        async with async_session_maker() as session:
            async with session.begin():
                query = select(cls.model).where(cls.model.pubchem_id == pubchem_id)
                molecule = await session.execute(query)

                return molecule.scalar_one_or_none()

    @classmethod
    async def delete_molecule_by_pubchem_id(cls, pubchem_id: str):
        async with async_session_maker() as session:
            async with session.begin():
                query = select(cls.model).filter_by(pubchem_id=pubchem_id)
                result = await session.execute(query)
                molecule_to_delete = result.scalar_one_or_none()

                if molecule_to_delete is None:
                    return None

                await session.execute(delete(cls.model).where(cls.model.pubchem_id == pubchem_id))

                await session.commit()
                return pubchem_id
