# TODO: Create plotly visualization of Venn Diagram of the read in data

import streamlit as st

def main():
    st.set_page_config(
        page_title="Spectrum Viewer",
        layout="wide",
    )
    st.sidebar.title("Upload files")
    uploaded_files = st.sidebar.file_uploader(
        "Choose tsv files", accept_multiple_files=True, type="tsv", key="uploader"
    )

# Run the main function
if __name__ == "__main__":
    main()
