"""Модуль для взаимодействия с пользователем"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from protein import (
    get_results_from_ncbi,
    multiple_alignment,
    plot_multiple_alignment,
    filter_sequences,
)


def main():
    """
    Основная функция проекта.

    1. Получение запроса от пользователя
    2. Получение ответа от NCBI
    3. Фильтрация выбросов по длине последовательности
    4. Сохранение данных в .csv
    5. Вывод изображения множественного выравнивания и расчёт матрицы попарных расстояний
    6. Сохранение матрицы попарных расстояний в .png и в .csv
    """
    # get user input
    while True:
        query = input("Enter your NCBI database search query:\n")
        if len(query) == 0:
            print("Error: empty query. Try again")
        else:
            break

    # get results from NCBI database
    max_records_return = 15
    results = get_results_from_ncbi(query, max_records_return)
    if len(results) == 0:
        print("Error: no results found")
        return

    # filter results by 1.5 IQR rule
    print(f"Length before filtering: {len(results)}")
    results = filter_sequences(results)
    print(f"Length after filtering: {len(results)}")

    # save info as .csv
    protein_df = pd.DataFrame(results)
    print(f"{query} search results:\n", protein_df)
    query_str = "_".join([x for x in query.split()])
    protein_df.to_csv(f"{query_str}_search_info.csv")

    # get multiple alignment data and plot
    alignment, distance_matrix = multiple_alignment(results)

    plot_multiple_alignment(alignment)

    # save distance matrix as numbers
    np.savetxt(f"{query_str}_distance_matrix.csv", distance_matrix, delimiter=",")

    plt.imshow(distance_matrix)
    plt.title("Distance matrix")
    plt.colorbar()
    plt.savefig(f"{query_str}_distance_matrix.png", dpi=300, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    # if this code was run as 'python ./main.py'
    main()
