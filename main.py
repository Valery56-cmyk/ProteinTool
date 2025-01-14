"""Модуль для взаимодействия с пользователем"""

from Bio import Entrez
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config import ALIGNMENT_PATH, DIST_MATRIX_PATH, MAX_SEQUENCES_NCBI, USER_EMAIL
from database import Database
from protein import (
    get_results_from_ncbi_and_database,
    multiple_alignment,
    plot_multiple_alignment,
    filter_sequences,
    alignment_to_fasta
)
from utils import logger


def main():
    """
    Основная функция проекта.

    1. Получение запроса от пользователя
    2. Загрузка локальной базы данных
    3. Получение ответа от NCBI и базы данных
    4. Обновление локальной базы данных
    5. Фильтрация выбросов по длине последовательности
    6. Сохранение данных в .csv
    7. Вывод изображения множественного выравнивания и расчёт матрицы попарных расстояний
    8. Сохранение матрицы попарных расстояний в .png и в .csv
    """
    # Configure Entrez
    Entrez.email = USER_EMAIL

    # get user input
    while True:
        query = input("Enter your NCBI database search query:\n")
        if len(query) == 0:
            logger.error("Error: empty query. Try again")
        else:
            break

    # load local database
    database = Database()

    # get results from NCBI database and local database
    results = get_results_from_ncbi_and_database(query, MAX_SEQUENCES_NCBI, database=database)
    if len(results) == 0:
        logger.error("Error: no results found")
        return
    
    # filter results by 1.5 IQR rule
    logger.info(f"Length before filtering: {len(results)}")
    results = filter_sequences(results)
    logger.info(f"Length after filtering: {len(results)}")
    
    # update database
    database.save_protein_sequences(results)

    # save info as .csv
    protein_df = pd.DataFrame(results)
    logger.info(f"{query} search results:\n{protein_df}")
    query_str = "_".join([x for x in query.split()])
    protein_df.to_csv(f"{query_str}_search_info.csv")

    # get multiple alignment data and plot
    alignment, distance_matrix = multiple_alignment(results)
    alignment.ids = [x["id"] for x in results]
    alignment_to_fasta(alignment, ALIGNMENT_PATH)

    # save distance matrix as numbers
    np.savetxt(DIST_MATRIX_PATH, distance_matrix, delimiter=",")


if __name__ == "__main__":
    # if this code was run as 'python ./main.py'
    main()
