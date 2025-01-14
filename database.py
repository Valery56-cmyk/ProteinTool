"""Модуль с классом базы данных белковых последовательностей"""

from contextlib import contextmanager
from datetime import datetime
from typing import List, Dict

from sqlalchemy import create_engine
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm import sessionmaker

from db_model import Base, ProteinSequence
from config import DB_PATH
from utils import logger


class Database:
    """
    Класс базы данных белковых последовательностей
    """

    def __init__(self):
        self.engine = create_engine(f"sqlite:///{DB_PATH}")
        Base.metadata.create_all(self.engine)
        self.SessionLocal = sessionmaker(bind=self.engine)

    @contextmanager
    def get_session(self):
        """
        Возвращает объект context manager для атомарной работы с базой данных.

        Атомарность заключается в том, что при ошибке при какой-либо операции, данная операция целиком отменяется (session.rollback),
        то есть данные в базе всегда валидны.
        """
        session = self.SessionLocal()
        try:
            yield session
            session.commit()
        except Exception as e:
            session.rollback()
            raise e
        finally:
            session.close()

    def save_protein_sequences(self, protein_sequences: List[Dict]):
        """
        Добавляет список белковых последовательностей в базу данных
        """
        with self.get_session() as session:
            for protein_sequence_data in protein_sequences:
                try:
                    # Create protein sequence instance
                    protein_sequence = self._protein_sequence_from_dict(
                        protein_sequence_data
                    )

                    # Merge protein sequence (update if exists, insert if new)
                    session.merge(protein_sequence)

                except SQLAlchemyError as e:
                    logger.error(
                        f"Error saving protein sequence {protein_sequence_data['id']}: {str(e)}"
                    )
                    raise

    def get_protein_sequences(self, limit: int = None) -> List[Dict]:
        """
        Возвращает `limit` первых записей о белковых последовательностях из базы данных
        """
        with self.get_session() as session:
            query = session.query(ProteinSequence)
            if limit:
                query = query.limit(limit)

            protein_sequences = query.all()
            return [
                self._protein_sequence_to_dict(protein_sequence)
                for protein_sequence in protein_sequences
            ]

    def get_protein_sequence_by_id(self, protein_sequence_id: str) -> Dict:
        """
        Возвращает белковую последовательность, если `protein_sequence_id` есть в базе.
        Иначе возвращает `None`
        """
        with self.get_session() as session:
            protein_sequence = (
                session.query(ProteinSequence)
                .filter(ProteinSequence.id == protein_sequence_id)
                .first()
            )
            if protein_sequence:
                return self._protein_sequence_to_dict(protein_sequence)
            return None

    def get_protein_sequence_by_organism(self, protein_sequence_organism: str) -> Dict:
        """
        Возвращает белковые последовательности, если `protein_sequence_organism` есть в базе.
        Иначе возвращает `None`
        """
        with self.get_session() as session:
            protein_sequences = (
                session.query(ProteinSequence).filter(
                    ProteinSequence.organism == protein_sequence_organism
                )
            ).all()
            if protein_sequences:
                return [self._protein_sequence_to_dict(x) for x in protein_sequences]
            return None

    @staticmethod
    def _protein_sequence_to_dict(protein_sequence: ProteinSequence) -> Dict:
        """
        Конвертирует объект `protein_sequence` типа `ProteinSequence` в словарь
        """
        return {
            "id": protein_sequence.id,
            "sequence": protein_sequence.sequence,
            "organism": protein_sequence.organism,
            "description": protein_sequence.description,
            "last_access_date": protein_sequence.last_access_date.strftime("%Y-%m-%d"),
        }

    @staticmethod
    def _protein_sequence_from_dict(protein_sequence_data: Dict) -> ProteinSequence:
        """
        Конвертирует словарь `protein_sequence_data` в объект типа `ProteinSequence`
        """
        return ProteinSequence(
            id=protein_sequence_data["id"],
            sequence=protein_sequence_data["sequence"],
            organism=protein_sequence_data["organism"],
            description=protein_sequence_data.get("description"),
            last_access_date=datetime.strptime(
                protein_sequence_data.get(
                    "last_access_date", datetime.now().strftime("%Y-%m-%d")
                ),
                "%Y-%m-%d",
            ).date(),
        )

    def delete_protein_sequence(self, protein_sequence_id: str):
        """
        Удаляет белковую последовательность с `protein_sequence_id` из базы данных
        """
        with self.get_session() as session:
            protein_sequence = (
                session.query(ProteinSequence)
                .filter(ProteinSequence.id == protein_sequence_id)
                .first()
            )
            if protein_sequence:
                session.delete(protein_sequence)
