import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import os
from adjustText import adjust_text
from sys import platform
import numpy as np
from sklearn.linear_model import LinearRegression

# Determine platform for path
if "win" in platform:
    slash = "\\"
else:
    slash = "/"
# Fix import and download paths
global DATA_PATH
global DOWNLOAD_PATH
global ION_LABELS
DATA_PATH = os.path.join(os.path.dirname(__file__).rstrip("pages"), "data")
DOWNLOAD_PATH = os.path.join(os.path.dirname(__file__).rstrip("pages"), "plots")
# Determine ion labels of interest
ION_LABELS = ["b", "y", "[M"]


# Read in uploaded file
@st.cache
def read_tsv(file):
    df = pd.read_csv(f"{DATA_PATH}/{file}", sep="\t")
    return df


def get_peaks_of_interest(all_peaks, matching_ions, peptide_of_interest):
    # Get peaks of the peptide of interest
    run_matching_ions = matching_ions.loc[
        matching_ions["Peptide"] == peptide_of_interest
    ]
    run_all_peaks = all_peaks.loc[all_peaks["Peptide"] == peptide_of_interest]

    return [run_all_peaks, run_matching_ions]


def get_colors(target_ions, exp_matching_ions):
    red = (1, 0, 0, 1)  # b-ions
    green = (0.2, 1, 0.5, 1)  # Precursor
    blue = (0, 0, 1, 1)  # y-ions
    grey = (0.5, 0.5, 0.5, 0.3)  # Other ions
    colors = []
    for index, ion in target_ions.iterrows():
        if ion["Is_matching_ion"]:
            label = exp_matching_ions["Ion_name"][
                exp_matching_ions["Ion_ID"] == index
            ].values[0]
            color = (
                red
                if "b" in label
                else blue
                if "y" in label
                else green
                if "[M" in label
                else grey
            )
        else:
            color = grey

        colors.append(color)

    return colors


# Compute relative intensity of each ion
def get_relative_intensities(intensities):
    return [intensity / max(intensities) for intensity in intensities]


# Compute correlation between two spectra
def spectrum_correlation(matching_ions_target1, matching_ions_target2):
    """
    Fit linear regression model onto overlapping matching ions between two spectra.
    Process two matplotlib ax.text objects and return the goodness of fit.
    :param matching_ions_target1:
    :param matching_ions_target2:
    :return: r2
    """
    # Strip ion mass and get overlapping ion labels
    overlapping_ions = set([ion.get_text().split('\n')[0] for ion in matching_ions_target1]) & \
                       set([ion.get_text().split('\n')[0] for ion in matching_ions_target2])
    if len(overlapping_ions) == 0:
        st.warning("No overlapping ions found. Cannot compute correlation.")
        return 0

    # Access labels and intensities of matplotlib text objects
    target1_intensities = np.array([ion.get_position()[1] for ion in matching_ions_target1 if ion.get_text().split('\n')[0] in overlapping_ions])
    target2_intensities = np.array([ion.get_position()[1] for ion in matching_ions_target2 if ion.get_text().split('\n')[0] in overlapping_ions])
    # Fit a linear regression model
    reg = LinearRegression().fit(target1_intensities.reshape(-1, 1), target2_intensities)

    # predict y values using the model ->  For future visualization
    # y_pred = reg.predict(x.reshape(-1, 1))

    # calculate R-squared (coefficient of determination) as a measure of goodness of fit
    r2 = reg.score(target1_intensities.reshape(-1, 1), target2_intensities)

    return r2


# Plot single spectrum
def plot_spectrum(df_all_peaks, df_matching_ions, peptide_of_interest, column):
    # Get peaks of interest
    df_target_all_ions, df_target_matching_ions = get_peaks_of_interest(
        df_all_peaks, df_matching_ions, peptide_of_interest
    )

    for replicate in df_target_all_ions["filename"].unique():
        rep_target_all_ions = df_target_all_ions[
            df_target_all_ions["filename"] == replicate
        ]
        # Compute relative intensities
        rep_target_all_ions["Relative_Intensity"] = get_relative_intensities(
            rep_target_all_ions["Intensity"]
        )
        # Get colors for each annotated ion. Blue for b-ions, red for y-ions, green for precursor
        rep_target_all_ions["color"] = get_colors(
            rep_target_all_ions, df_target_matching_ions
        )

        # Plot spectrum
        fig, ax = plt.subplots(figsize=(15, 10))
        # Plot experimental spectrum
        ax.bar(
            rep_target_all_ions["Experimental_mass"],
            rep_target_all_ions["Relative_Intensity"],
            color=rep_target_all_ions["color"],
            width=st.session_state["bar_width_slider"],
        )
        # Get ion labels
        ion_labels = []
        for i, ion in rep_target_all_ions.iterrows():
            if ion["Is_matching_ion"]:
                x = ion["Experimental_mass"]
                y = ion["Relative_Intensity"]
                label = df_target_matching_ions["Ion_name"][
                    df_target_matching_ions["Ion_ID"] == i
                ].values[0]
                if (
                    any(valid_label in label for valid_label in ION_LABELS)
                    and y >= st.session_state["show_label_slider"]
                ):
                    label_and_mass = f"{label}\n{x:.3f}"
                    ion_labels.append(
                        ax.text(
                            x,
                            y,
                            label_and_mass,
                            color=ion["color"],
                            fontsize=st.session_state["tick_label_size_slider"],
                        )
                    )
        # Adjust text to avoid overlapping
        adjust_text(
            ion_labels,
            only_move={"points": "y", "texts": "y"},
            arrowprops=dict(arrowstyle="-", color="black", lw=1),
            color="red",
        )

        min_x = round(min(rep_target_all_ions["Experimental_mass"]), 0) - 10
        max_x = round(max(rep_target_all_ions["Experimental_mass"]), 0) + 10
        ax.set_xlim([min_x, max_x])
        ax.set_ylabel(
            "Relative Intensity", fontsize=st.session_state["legend_label_size_slider"]
        )
        ax.set_xlabel("m/z", fontsize=st.session_state["legend_label_size_slider"])
        ax.tick_params(labelsize=st.session_state["tick_label_size_slider"])
        ax.tick_params(labelsize=st.session_state["tick_label_size_slider"])
        ax.title.set_text(f"MS2 Spectrum of {peptide_of_interest} in {replicate}")

        # plot this to streamlit
        with column:
            st.pyplot(fig)
            st.button(
                "Download",
                key=f"download_button_peptide_{peptide_of_interest}_{replicate}",
                on_click=handle_download_button,
                args=(fig, peptide_of_interest, replicate),
            )


# Plot mirror plot
def plot_mirror_plot(
    df1_all_peaks,
    df1_matching_ions,
    df2_all_peaks,
    df2_matching_ions,
    peptide_of_interest,
    column,
):
    # Get peaks of interest
    df_target1_all_ions, df_target1_matching_ions = get_peaks_of_interest(
        df1_all_peaks, df1_matching_ions, peptide_of_interest
    )
    df_target2_all_ions, df_target2_matching_ions = get_peaks_of_interest(
        df2_all_peaks, df2_matching_ions, peptide_of_interest
    )

    for i, replicate_target1 in enumerate(df_target1_all_ions["filename"].unique()):
        for y, replicate_target2 in enumerate(df_target2_all_ions["filename"].unique()):
            rep_target1_all_ions = df_target1_all_ions[
                df_target1_all_ions["filename"] == replicate_target1
            ]
            rep_target2_all_ions = df_target2_all_ions[
                df_target2_all_ions["filename"] == replicate_target2
            ]
            # Compute relative intensities
            rep_target1_all_ions["Relative_Intensity"] = get_relative_intensities(
                rep_target1_all_ions["Intensity"]
            )
            rep_target2_all_ions["Relative_Intensity"] = get_relative_intensities(
                rep_target2_all_ions["Intensity"]
            )
            # Get colors for each annotated ion. Blue for b-ions, red for y-ions, green for precursor
            rep_target1_all_ions["color"] = get_colors(
                rep_target1_all_ions, df_target1_matching_ions
            )
            rep_target2_all_ions["color"] = get_colors(
                rep_target2_all_ions, df_target2_matching_ions
            )

            # Plot mirror plot
            fig, ax = plt.subplots(figsize=(15, 10))
            # Plot experimental spectrum
            ax.bar(
                rep_target1_all_ions["Experimental_mass"],
                rep_target1_all_ions["Relative_Intensity"],
                color=rep_target1_all_ions["color"],
                width=st.session_state["bar_width_slider"],
            )
            ax.bar(
                rep_target2_all_ions["Experimental_mass"],
                rep_target2_all_ions["Relative_Intensity"] * -1,
                color=rep_target2_all_ions["color"],
                width=st.session_state["bar_width_slider"],
            )
            # Add horizontal line at 0
            ax.axhline(0, color="grey", linewidth=0.8, alpha=0.3)
            # Get labels
            ion_labels_target1 = []
            ion_labels_target2 = []
            for i, ion in rep_target1_all_ions.iterrows():
                if ion["Is_matching_ion"]:
                    x = ion["Experimental_mass"]
                    y = ion["Relative_Intensity"]
                    label = df_target1_matching_ions["Ion_name"][
                        df_target1_matching_ions["Ion_ID"] == i
                    ].values[0]
                    if (
                        any(valid_label in label for valid_label in ION_LABELS)
                        and y >= st.session_state["show_label_slider"]
                    ):
                        label_and_mass = f"{label}\n{x:.3f}"
                        ion_labels_target1.append(
                            ax.text(
                                x,
                                y,
                                label_and_mass,
                                color=ion["color"],
                                fontsize=st.session_state["tick_label_size_slider"],
                            )
                        )

            for j, ion in rep_target2_all_ions.iterrows():
                if ion["Is_matching_ion"]:
                    x = ion["Experimental_mass"]
                    y = ion["Relative_Intensity"] * -1
                    label = df_target2_matching_ions["Ion_name"][
                        df_target2_matching_ions["Ion_ID"] == j
                    ].values[0]
                    if (
                        any(valid_label in label for valid_label in ION_LABELS)
                        and y <= st.session_state["show_label_slider"] * -1
                    ):
                        label_and_mass = f"{label}\n{x:.3f}"
                        ion_labels_target2.append(
                            ax.text(
                                x,
                                y,
                                label_and_mass,
                                color=ion["color"],
                                fontsize=st.session_state["tick_label_size_slider"],
                            )
                        )
            # Compute correlation coefficient of matching ions
            #cor_coeff =
            ion_labels = ion_labels_target1 + ion_labels_target2
            adjust_text(
                ion_labels,
                only_move={"points": "y", "texts": "y"},
                arrowprops=dict(arrowstyle="-", color="black", lw=1),
                color="red",
            )

            min_x = (
                round(
                    min(
                        min(rep_target1_all_ions["Experimental_mass"]),
                        min(df_target2_all_ions["Experimental_mass"]),
                    ),
                    0,
                )
                - 10
            )
            max_x = (
                round(
                    max(
                        max(rep_target1_all_ions["Experimental_mass"]),
                        max(df_target2_all_ions["Experimental_mass"]),
                    ),
                    0,
                )
                + 10
            )
            ax.title.set_text(
                f"Mirror plot of {peptide_of_interest} in {replicate_target1} and {replicate_target2}"
            )
            ax.set_xlim([min_x, max_x])
            ax.set_ylabel(
                "Relative Intensity",
                fontsize=st.session_state["legend_label_size_slider"],
            )
            ax.set_xlabel("m/z", fontsize=st.session_state["legend_label_size_slider"])
            ax.tick_params(labelsize=st.session_state["tick_label_size_slider"])
            ax.tick_params(labelsize=st.session_state["tick_label_size_slider"])

            # Add goodness of fit to lower right of plot
            ax.text(0.95, 0.05, f"R² = {spectrum_correlation(ion_labels_target1, ion_labels_target2):.2f}",
                    transform=ax.transAxes,
                    fontsize=st.session_state["tick_label_size_slider"],
                    verticalalignment='bottom',
                    horizontalalignment='right')

            # plot this to streamlit
            with column:
                st.pyplot(fig)
                # Add download button
                st.button(
                    "Download",
                    key=f"download_button_peptide_{peptide_of_interest}_{replicate_target1}_{replicate_target2}",
                    on_click=handle_download_button,
                    args=(
                        fig,
                        peptide_of_interest,
                        replicate_target1,
                        replicate_target2,
                    ),
                )


def handle_download_button(fig, peptide_of_interest, target1, target2=None):
    if target2 is not None:
        fig.savefig(
            f"{DOWNLOAD_PATH}/{peptide_of_interest}_{target1}_{target2}.{st.session_state['export_format']}",
            dpi=300,
        )
        st.success(
            f"Downloaded to {DOWNLOAD_PATH}/{peptide_of_interest}_{target1}_{target2}.{st.session_state['export_format']}"
        )
    else:
        fig.savefig(
            f"{DOWNLOAD_PATH}/{peptide_of_interest}_{target1}.{st.session_state['export_format']}",
            dpi=300,
        )
        st.success(
            f"Downloaded to {DOWNLOAD_PATH}/{peptide_of_interest}_{target1}.{st.session_state['export_format']}"
        )


def main():
    # TODO: set page icon
    st.set_page_config(
        page_title="Spectrum Viewer",
        layout="wide",
        page_icon="simple-icons:soundcharts",
    )
    st.sidebar.title("Upload files")
    uploaded_files = st.sidebar.file_uploader(
        "Choose tsv files", accept_multiple_files=True, type="tsv", key="uploader"
    )

    col_plot, col_settings = st.columns([3, 1])
    with col_settings:
        st.title("Settings")
        st.slider(
            "Threshold for displaying ion annotation",
            0.0,
            1.0,
            value=0.1,
            step=0.05,
            key="show_label_slider",
        )
        st.slider(
            "Axis label size", 0, 30, value=15, step=1, key="legend_label_size_slider"
        )
        st.slider(
            "Tick label size", 0, 20, value=10, step=1, key="tick_label_size_slider"
        )
        st.slider("Bar width", 0.0, 4.0, value=2.0, step=0.1, key="bar_width_slider")
        st.selectbox(
            "Export format", ["svg", "png", "jpeg", "pdf"], key="export_format"
        )
    with col_plot:
        st.title("Plot")
        # Only start if files are uploaded
        if len(st.session_state["uploader"]) != 0:
            uploaded_file_names = [file.name for file in uploaded_files]

            # Select files of Spectrum 1
            st.sidebar.text("Select files of Spectrum 1")
            st.sidebar.selectbox(
                "All peaks", uploaded_file_names, key="file1_all_peaks"
            )
            st.sidebar.selectbox(
                "Matching Ions", uploaded_file_names, key="file1_matching_ions"
            )
            # Read files for spectrum 1
            df1_all_peaks = read_tsv(st.session_state["file1_all_peaks"])
            df1_matching_ions = read_tsv(st.session_state["file1_matching_ions"])

            # Select files of Spectrum 2 if checkbox is checked
            select_mirror_plot = st.sidebar.checkbox(
                "Compare with Spectrum 2", key="checkbox_mirror_plot"
            )

            if select_mirror_plot:
                st.sidebar.text("Select files of Spectrum 2")
                st.sidebar.selectbox(
                    "All peaks", uploaded_file_names, key="file2_all_peaks"
                )
                st.sidebar.selectbox(
                    "Matching ions", uploaded_file_names, key="file2_matching_ions"
                )
                # Read files for spectrum 2
                df2_all_peaks = read_tsv(st.session_state["file2_all_peaks"])
                df2_matching_ions = read_tsv(st.session_state["file2_matching_ions"])

                peptide_selection = set(df1_all_peaks["Peptide"]) & set(df2_all_peaks["Peptide"])

                # Select peptide of interest
                st.selectbox(
                    f"Choose a peptide out of {len(peptide_selection)}",
                    peptide_selection,
                    key="peptide_of_interest",
                )
                # Create a button that plots the mirror plot
                st.button(
                    "Plot",
                    key="button_mirror_plot",
                    on_click=plot_mirror_plot,
                    args=(
                        df1_all_peaks,
                        df1_matching_ions,
                        df2_all_peaks,
                        df2_matching_ions,
                        st.session_state["peptide_of_interest"],
                        col_plot,
                    ),
                )
            else:
                peptide_selection = set(df1_all_peaks["Peptide"])
                # Select peptide of interest
                st.selectbox(
                    f"Choose a peptide out of {len(peptide_selection)}",
                    peptide_selection,
                    key="peptide_of_interest",
                )
                # Create a button that plots the spectrum
                st.button(
                    "Plot",
                    key="button_plot_spectrum",
                    on_click=plot_spectrum,
                    args=(
                        df1_all_peaks,
                        df1_matching_ions,
                        st.session_state["peptide_of_interest"],
                        col_plot,
                    ),
                )

    # Create a button that clears out the right pane
    if st.sidebar.button("Clear"):
        st.empty()


# Run the main function
if __name__ == "__main__":
    main()
