"""Модуль с описанием структуры базы данных белковых последовательностей"""

from sqlalchemy import Column, String, Date
from sqlalchemy.orm import declarative_base

Base = declarative_base()


class ProteinSequence(Base):
    """
    Класс белковой последовательности - одна строка в базе данных

    `id` - уникальный идентификатор записи

    `sequence` - последовательность аминокислот

    `organism` - организм

    `description` - описание, дополнительная информация о записи

    `last_access_date` - дата последнего доступа
    """

    __tablename__ = "protein_sequences"

    id = Column(String, unique=True, index=True, primary_key=True)
    sequence = Column(String, nullable=False)
    organism = Column(String, index=True)
    description = Column(String)
    last_access_date = Column(Date, nullable=False)
