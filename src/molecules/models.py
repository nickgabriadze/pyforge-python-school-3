from sqlalchemy.orm import Mapped
from ..db.database import Base, str_uniq, pubchem_pk


class Molecule(Base):
    pubchem_id: Mapped[pubchem_pk]
    smiles: Mapped[str_uniq]

    def __str__(self):
        return (
            f"{self.__class__.__name__}(id={self.pubchem_id}, "
            f"smiles={self.smiles!r})"
        )

    def __repr__(self):
        return str(self)

    def getSmiles(self):
        return self.smiles

    def getPubchemID(self):
        return self.pubchem_id
