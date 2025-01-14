"""Модуль для визуализации данных в виде дашборда"""

import dash
from dash import dash, html
from dash import html, dcc, callback, Output, Input
import dash_bio as dashbio
import numpy as np
import pandas as pd
import plotly.express as px

from database import Database
from config import ALIGNMENT_PATH, DIST_MATRIX_PATH


class Dashboard:
    """Класс дашборда"""

    def __init__(self):
        """Инициализация данных для визуализации в виде pandas-таблицы"""
        Dashboard._database = Database()
        Dashboard.df = pd.DataFrame(
            Dashboard._database.get_protein_sequences(limit=None),
            columns=["id", "sequence", "organism", "description", "last_access_date"],
        )
        Dashboard.matrix = np.genfromtxt(DIST_MATRIX_PATH, delimiter=",")

        with open(ALIGNMENT_PATH, "r") as f:
            Dashboard.alignment_data = f.readlines()

        self.fig_organism = self.get_organism_hist_fig()
        self.fig_sequence = self.get_sequence_length_distr_fig()
        self.fig_last_access = self.get_last_access_date_distr_fig()
        self.fig_matrix = self.get_dist_matrix_heatmap_fig(matrix=self.matrix)

    def get_organism_hist_fig(self):
        """Визуализация категориальной гистограммы распределения последовательностей по организмам"""
        # get data
        organism_counts = self.df["organism"].value_counts().reset_index()
        organism_counts.columns = ["organism", "count"]

        # get histogram
        fig_organism = px.histogram(
            organism_counts,
            x="organism",
            y="count",
            title="Organism Distribution",
            labels={"organism": "Organism", "count": "Number of Sequences"},
            hover_data=["organism", "count"],
        )
        fig_organism.update_layout(barmode="stack")
        return fig_organism

    def get_sequence_length_distr_fig(self, bin_size: int = 10):
        """Визуализация гистограммы распределения последовательностей по длине последовательности"""
        # get data
        self.df["sequence_length"] = self.df["sequence"].apply(len)

        # get histogram
        fig_sequence = px.histogram(
            self.df,
            x="sequence_length",
            title="Distribution of Sequence Lengths",
            labels={"sequence_length": "Sequence Length"},
            hover_data=["sequence_length"],
        )
        fig_sequence.update_traces(xbins_size=bin_size)
        return fig_sequence

    def get_last_access_date_distr_fig(self):
        """Визуализация гистограммы распределения последовательностей по дате последнего обращения"""
        # get data
        self.df["last_access_date"] = pd.to_datetime(self.df["last_access_date"])

        # get histogram
        fig_last_access = px.histogram(
            self.df,
            x="last_access_date",
            title="Distribution of Last Access Dates",
            labels={"last_access_date": "Last Access Date"},
            hover_data=["last_access_date"],
        )
        fig_last_access.update_xaxes(rangeslider_visible=True)
        return fig_last_access

    def get_dist_matrix_heatmap_fig(self, matrix: np.ndarray = np.zeros((100, 100))):
        """Визуализация тепловой карты распределения расстояния между последовательностями"""
        fig_matrix = px.imshow(
            matrix,
            labels=dict(
                x="Protein sequence 1", y="Protein sequence 2", color="Distance"
            ),
            title="Distance Matrix Heatmap",
            color_continuous_scale="RdBu_r",
        )
        fig_matrix.update_layout(title_x=0.5)
        fig_matrix.update_layout(width=1000, height=1000)
        return fig_matrix

    def run_visualization(self, debug: bool = True):
        """Запуск интерактивной dash-визуализации"""
        organism_options = [{"label": "All", "value": "All"}] + [
            {"label": org, "value": org} for org in self.df["organism"].unique()
        ]

        app = dash.Dash(__name__)
        # all components: headers, plots, etc
        app.layout = html.Div(
            [
                html.H1("Protein Sequence Dashboard"),
                dcc.Dropdown(
                    id="organism-dropdown",
                    options=organism_options,
                    value="All",
                    clearable=False,
                ),
                dcc.Graph(id="organism-histogram", figure=self.fig_organism),
                dcc.Graph(id="sequence-length-histogram", figure=self.fig_sequence),
                dcc.Graph(id="last-access-date-histogram", figure=self.fig_last_access),
                #########################################
                html.H1("Last query results"),
                dashbio.AlignmentChart(
                    id="my-default-alignment-viewer",
                    data=self.alignment_data,
                    height=1200,
                    tilewidth=30,
                ),
                html.Div(id="default-alignment-viewer-output"),
                dcc.Graph(id="matrix-heatmap", figure=self.fig_matrix),
                # # show consensus sequence?
                # dashbio.SequenceViewer(
                #     id="default-sequence-viewer",
                #     sequence="MPASSQTAASERRNERVLTHEPPKLKTTLKTPPKQATNPLQFVKVGPCSLYRTAQEQLQKVQEVKKIKQEVRDDPEDWQSNLDNWKSSRRKRQEHIIERVVEVKKLELEEHDRQRRRNKTFSEMMEERGNRGRKLSISLAMYNDEDANDLSDLGIGTSSGKSSVSGDTHDDTHSVLSDRDSEIEKSHSDVDNAPDMTSSAATLTTTGTTVTTTSATATIATLTTTTTNISTMVVSTTAPTTTTGTTSTARKMFSGGFHSNLNHNQEYDSGTTATTSSPEPEEYTYEGAIRGYVSRVSQNIPRRSLTGVDSKLETAKSTINGSKTNVNDDGSKSPLSVVKVDILKRREVFEKASQKSNDNKANNRLSGDFTGTKSIKERLSSLERQKYETENSDKATNKTLNRLSGDMSSIRERLTHLEKQASERESKSSVHRKLSTEDLETGRPLRERLSTLEKYSSSDESSVPITTMESHSHNGDITARTIKDRLSALDPARGKETTDKRGPGHVGKHPLCFRDQENRVDISTPSERSSSPDSEYRVPRAVFHRSLDSLDADASSGPDTFERVQSLEELDYGRQYPASSSSAELLNDTDREDSGIHTADVSCSVSQADEPIDEEIVHHPAGVVVERREVTVEERVKATSIIQEDASTSNTTVEAPSDAVVATHSQFASHRETVAVNKERIAVEKADQEVEAAATACAEDPAAPAVESKQQEATETSTSVEMPAAAKKSDVEVDVHDVPNATVDAPTAKHPLERPTNLSLAKQEVGLADALNIPSPASPISKQILPLTLANDEEIIAGQPFLLNPPTSVEPPKEKPPPPPVDVSDDENPPPEPLKRLNSTRRIKKELRTRRSDFLGIEGINDDDLEPELTLTKPPDMAAILAEERRIEQLHRRSYDTDSNYEQDSSHERDSGVELGHVEDWAKQPVSPDMSQHSRQSSEPFGASVTSSEEDEITKKEREIIEVLEKEEQWRYGDNREYNRYCNHKYYNSELGERLAHKLRELEEEKMQLEREAILRCNEDAFRKRDDNSRKQEDASSHEQSQQPHDETAKQREIQARTTEEEAKYVEDKRKDEELETQSEQLRMQNEIVEKERRVADACRKEEMRIRTMESQIREQESKALVGAGLSSNEDEFPSGEVLRVERELLQLEQEELKRQRNNLAYREQKQLRLAEQLQEQWKSLQDVAQNSIKNTQQYKYHTPAVNYRSSMPDLQFQDVPRRRPPPPSIPLTKPRMLDQRQRDVTIRNSRIPSADSIPQQVDASIRESATATLGNHPGQQMSRQTLQALSAVPRPRIVQGDQWVQRRKSDVPRGAHDVNYQHWLIQEAEQRRINERNQRSPARKSQPHVTGTSVPYNAPIRTDSKPLPDSIIQTLTQRVQNKTQEKPLSTRRRPEQVYNQEQHLTVQHQSPQYTLQSQKTQPAPIANNNGNQEKMLSVSGKKKCSHCGEELGRGAAMIIESLRLFYHMECFKCCVCHVRLGDGLIGTDVRVRNHKLHCHNCYSSDDGVKFSCV",
                #     title="Example",
                #     badge=False,
                #     charsPerLine=100,
                # ),
                # html.Div(id="default-sequence-viewer-output"),
            ]
        )
        app.run_server(debug=debug)

    @callback(
        Output("organism-histogram", "figure"),
        Output("sequence-length-histogram", "figure"),
        Output("last-access-date-histogram", "figure"),
        Input("organism-dropdown", "value"),
    )
    def update_graphs_on_selected_organism(selected_organism):
        """Обновление визуализации при изменении выбранного организма"""
        if selected_organism == "All":
            filtered_df = Dashboard.df.copy()
        else:
            filtered_df = Dashboard.df[Dashboard.df["organism"] == selected_organism]

        # Organism histogram
        organism_counts = filtered_df["organism"].value_counts().reset_index()
        organism_counts.columns = ["organism", "count"]
        fig_organism = px.histogram(
            organism_counts,
            x="organism",
            y="count",
            title="Organism Distribution",
            labels={"organism": "Organism", "count": "Number of Sequences"},
            hover_data=["organism", "count"],
        )
        fig_organism.update_layout(barmode="stack")

        # Sequence length distribution
        fig_sequence = px.histogram(
            filtered_df,
            x="sequence_length",
            title="Distribution of Sequence Lengths",
            labels={"sequence_length": "Sequence Length"},
            hover_data=["sequence_length"],
        )
        fig_sequence.update_traces(xbins_size=10)

        # Last access date distribution
        fig_last_access = px.histogram(
            filtered_df,
            x="last_access_date",
            title="Distribution of Last Access Dates",
            labels={"last_access_date": "Last Access Date"},
            hover_data=["last_access_date"],
        )
        fig_last_access.update_xaxes(rangeslider_visible=True)

        return fig_organism, fig_sequence, fig_last_access


if __name__ == "__main__":
    dashboard = Dashboard()
    dashboard.run_visualization()
