"""Модуль для обработки белковых последовательностей"""

from biotite.sequence import ProteinSequence
import biotite.sequence.align as align
from Bio import Entrez, SeqIO
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from plot_alignment import plot_alignment_shapes


def get_results_from_ncbi(query: str, query_size: int = 10) -> dict:
    """
    Функция, которая отправляет запрос в базу белковых последовательностей NCBI.

    1. Entrez.esearch - поиск ID последовательностей по запросу пользователя.
    2. Entrez.efetch - поиск информации про последовательности из пункта 1.
    """
    result_sequences = []
    Entrez.email = "your.email@uc.cl"

    # get protein ids
    handle = Entrez.esearch(db="protein", retmax=query_size, term=query, idtype="acc")
    record = Entrez.read(handle)
    handle.close()

    # get protein data for selected ids
    handle = Entrez.efetch(db="protein", rettype="fasta", id=record["IdList"])
    for seq_record in SeqIO.parse(handle, "fasta"):
        result_sequences.append(seq_record.__dict__)
    handle.close()

    return result_sequences


def filter_sequences(protein_info_list: list) -> list:
    """
    Функция для фильтрации последовательностей по длине согласно 1.5 * IQR rule.

    Принимает список последовательностей, возвращает список отфильтрованных последовательностей.

    Статья с теорией:
    https://builtin.com/articles/1-5-iqr-rule
    """
    # quartiles - Q1, Q2, Q3
    # Q1 - 25%
    # Q2 - 50% (median)
    # Q3 - 75%
    # IQR = InterQuartile Range = Q3 - Q1
    # 1.5 * IQR rule:
    # [low, high]
    # low = Q1 - 1.5 * IQR
    # high = Q3 + 1.5 * IQR
    # 3 * sigma rule: [-3 * sigma, 3 * sigma] = 99.73% confidence

    Q3 = np.quantile([len(x["_seq"]) for x in protein_info_list], q=0.75)
    Q1 = np.quantile([len(x["_seq"]) for x in protein_info_list], q=0.25)

    IQR = Q3 - Q1
    low = Q1 - 1.5 * IQR
    high = Q3 + 1.5 * IQR

    filtered_sequences = []
    for protein in protein_info_list:
        if low <= len(protein["_seq"]) <= high:
            filtered_sequences.append(protein)
    return filtered_sequences


def multiple_alignment(protein_info_list: list) -> tuple:
    """
    Функция для множественного выравнивания последовательностей.

    Возвращает данные выравнивания и матрицу попарных расстояний.
    """
    # get result sequences
    seqs = []
    for protein in protein_info_list:
        seqs.append(ProteinSequence(protein["_seq"]))

    # Perform a multiple sequence alignment
    matrix = align.SubstitutionMatrix.std_protein_matrix()
    alignment, order, _, distance_matrix = align.align_multiple(seqs, matrix)
    # Order alignment according to guide tree
    alignment = alignment[:, order.tolist()]
    return alignment, distance_matrix


def plot_multiple_alignment(alignment):
    """Функция для построения графика множественного выравнивания последовательностей."""
    fig, ax = plt.subplots()

    plot_alignment_shapes(
        ax, alignment, labels=None, symbols_per_line=len(alignment), symbol_size=8
    )

    # ax.get_figure().patch.set_facecolor("#181818")

    ax.set_ylabel("Sequence", color="black")
    ax.set_title("Multiple sequence alignment", color="black")
    fig.tight_layout()

    # Adjust the bottom size according to the requirement of the user
    plt.subplots_adjust(bottom=0.25)

    # Choose the Slider color
    slider_color = "grey"

    min_x, max_x = ax.get_xlim()
    step = 10

    # Set the axis and slider position in the plot
    axis_position = plt.axes([0.2, 0.1, 0.65, 0.03], facecolor=slider_color)
    slider_position = Slider(axis_position, "X", min_x, max_x - step, valinit=min_x)
    ax.axis([min_x, step, None, None])

    # update() function to change the graph when the slider is in use
    def update(_):
        """Функция, которая вызывается при изменении положения ползунка"""
        pos = slider_position.val
        ax.axis([pos, pos + step, None, None])
        fig.canvas.draw_idle()

    # update function called using on_changed() function
    slider_position.on_changed(update)

    plt.show()
