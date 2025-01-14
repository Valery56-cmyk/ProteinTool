"""Модуль для обработки белковых последовательностей"""

from biotite.sequence import ProteinSequence
import biotite.sequence.align as align
from Bio import Entrez, SeqIO
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

from config import FASTA_PATH
from database import Database
from plot_alignment import plot_alignment_shapes
from utils import logger


def get_results_from_ncbi(
    query: str, query_size: int = 10, save_to_fasta: bool = False
) -> dict:
    """
    Функция, которая отправляет запрос в базу белковых последовательностей NCBI.

    1. Entrez.esearch - поиск ID последовательностей по запросу пользователя.
    2. Entrez.efetch - поиск информации про последовательности из пункта 1.
    """
    result_sequences = []

    # get protein ids
    handle = Entrez.esearch(db="protein", retmax=query_size, term=query, idtype="acc")
    record = Entrez.read(handle)
    handle.close()
    logger.info(f"Found {len(record['IdList'])} sequence IDs")

    # get protein data for selected ids
    handle = Entrez.efetch(db="protein", rettype="fasta", id=record["IdList"])

    # save to FASTA file if needed
    if save_to_fasta:
        sequences = SeqIO.parse(handle, "fasta")
        with open(FASTA_PATH, "w") as protein_output_handle:
            SeqIO.write(sequences, protein_output_handle, "fasta")

    for seq_record in SeqIO.parse(handle, "fasta"):
        result_sequences.append(seq_record.__dict__)
    handle.close()

    return result_sequences


def get_results_from_ncbi_and_database(
    query: str,
    query_size: int = 10,
    save_to_fasta: bool = False,
    database: Database = None,
) -> dict:
    """
    Функция, которая отправляет запрос в базу белковых последовательностей NCBI.

    1. Entrez.esearch - поиск ID последовательностей по запросу пользователя.
    2. Поиск ID последовательностей в базе данных.
    3. Entrez.efetch - поиск информации про оставшиеся последовательности.

    Также обновляет локальную базу данных.
    """
    result_sequences = []

    # get protein ids
    handle = Entrez.esearch(db="protein", retmax=query_size, term=query, idtype="acc")
    record = Entrez.read(handle)
    handle.close()

    ids_to_efetch = []
    if database:
        # if database usage is enabled
        for sequence_id in record["IdList"]:
            database_protein_sequence = database.get_protein_sequence_by_id(sequence_id)
            if database_protein_sequence:
                # found sequence in database => add to results
                result_sequences.append(database_protein_sequence)
                logger.info(f"Found {database_protein_sequence['id']} in database")
            else:
                # add sequence id for Entrez search query
                ids_to_efetch.append(sequence_id)
    else:
        # copy list of ids
        ids_to_efetch = record["IdList"][::]

    logger.info(f"Will search for {len(ids_to_efetch)} sequences in NCBI")

    # get protein data for selected ids
    handle = Entrez.efetch(db="protein", rettype="fasta", id=ids_to_efetch)

    # save to FASTA file if needed
    if save_to_fasta:
        sequences = SeqIO.parse(handle, "fasta")
        with open(FASTA_PATH, "w") as protein_output_handle:
            SeqIO.write(sequences, protein_output_handle, "fasta")

    for seq_record in SeqIO.parse(handle, "fasta"):
        logger.info(f"{seq_record = }")

        left = seq_record.description.find("[")
        right = seq_record.description.rfind("]")
        organism = ""
        if (left < right) and left != -1 and right != -1:
            organism = seq_record.description[left + 1 : right]
        result_sequences.append(
            {
                "id": seq_record.id,
                "sequence": str(seq_record._seq),
                "description": seq_record.description,
                "organism": organism,
            }
        )
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
    # q1 - 25%
    # q2 - 50% (median)
    # q3 - 75%
    # iqr = InterQuartile Range = Q3 - Q1
    # 1.5 * IQR rule:
    # [low, high]
    # low = Q1 - 1.5 * IQR
    # high = Q3 + 1.5 * IQR
    # 3 * sigma rule: [-3 * sigma, 3 * sigma] = 99.73% confidence

    Q3 = np.quantile([len(x["sequence"]) for x in protein_info_list], q=0.75)
    Q1 = np.quantile([len(x["sequence"]) for x in protein_info_list], q=0.25)

    IQR = Q3 - Q1
    low = Q1 - 1.5 * IQR
    high = Q3 + 1.5 * IQR

    filtered_sequences = []
    for protein in protein_info_list:
        if low <= len(protein["sequence"]) <= high:
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
        seqs.append(ProteinSequence(protein["sequence"]))

    # perform a multiple sequence alignment
    matrix = align.SubstitutionMatrix.std_protein_matrix()
    alignment, order, _, distance_matrix = align.align_multiple(seqs, matrix)
    # order alignment according to guide tree
    alignment = alignment[:, order.tolist()]
    return alignment, distance_matrix


def plot_multiple_alignment(alignment):
    """
    Функция для построения графика множественного выравнивания последовательностей.
    """
    fig, ax = plt.subplots()

    plot_alignment_shapes(
        ax, alignment, labels=None, symbols_per_line=len(alignment), symbol_size=8
    )

    # ax.get_figure().patch.set_facecolor("#181818")

    ax.set_ylabel("Sequence", color="black")
    ax.set_title("Multiple sequence alignment", color="black")
    fig.tight_layout()

    # adjust the bottom size according to the requirement of the user
    plt.subplots_adjust(bottom=0.25)

    # choose the Slider color
    slider_color = "grey"

    min_x, max_x = ax.get_xlim()
    step = 10

    # set the axis and slider position in the plot
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


def alignment_to_fasta(alignment, output_file, remove_gaps=False):
    """
    Конвертирует данные множественного выравнивания в FASTA файл
    """
    aligned_sequences = alignment.get_gapped_sequences()
    try:
        headers = alignment.ids
    except AttributeError:
        headers = [f"seq_{i+1}" for i in range(len(aligned_sequences))]

    with open(output_file, "w") as fasta_file:
        for header, seq in zip(headers, aligned_sequences):
            # Ensure the header starts with '>'
            if not header.startswith(">"):
                header = ">" + header
            # Convert sequence to string
            seq_str = str(seq)
            # Optionally remove gaps
            if remove_gaps:
                seq_str = seq_str.replace("-", "")
            # Write header and sequence
            fasta_file.write(f"{header}\n")
            fasta_file.write(f"{seq_str}\n")
