"""
This is the script the streamlit app uses
https://nanoparticleatomcounter.streamlit.app/
"""
import contextlib
import os
import tempfile
from pathlib import Path
import subprocess
import sys
import pandas as pd
import streamlit as st
from NanoparticleAtomCounter.cli.atom_count import main as atom_counter


st.set_page_config(page_title="Nanoparticle Atom Counter", page_icon="🧮")

with st.sidebar:
    st.header("Resources")

    SAMPLE_CSV = Path("docs/sample_input.csv").read_bytes()
    st.download_button(
        "Sample input (.csv)",
        SAMPLE_CSV,
        file_name="sample_input.csv",
        mime="text/csv",
    )

    st.image("docs/Acute_1.png", caption="θ < 90°", use_container_width=True)
    st.image("docs/Obtuse_1.png", caption="θ > 90°", use_container_width=True)
    st.image(
        "docs/Nanoparticle_Legend.png",
        caption="Definition of surface, interfacial, and perimeter atoms",
        use_container_width=True,
    )


st.title("Nanoparticle Atom Counter")

st.markdown(
    """
Upload a **.csv** containing the columns  
`r (A), R (A), Theta, Element, Interface Facet, Surface Facet`.

*Supply either **r** or **R** (if both are present, **r** is used).  
Interface Facet and Surface Facet are optional; leave blank if unknown.*

**Need a template or a visual guide?**  
A sample input file and explanatory diagrams are in the sidebar.
""",
    unsafe_allow_html=True,
)

mode = st.radio("**Select Counting Mode**", ("volume", "area"), horizontal=True)

file = st.file_uploader(
    "Drag-and-drop or browse for your input file",
    type=("csv", "xls", "xlsx"),
    accept_multiple_files=False,
)

if file is None:
    st.stop()  # wait for the user's input

if st.button("⚙️ Run calculation"):
    with st.spinner("Processing . . . "):

        # temp input file
        in_suffix = Path(file.name).suffix
        with tempfile.NamedTemporaryFile(delete=False, suffix=in_suffix) as tempin:
            tempin.write(file.getbuffer())
            tempin.flush()
            in_path = tempin.name

        # temp output file
        with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as tempout:
            out_path = tempout.name

        try:
            # run the calculation
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "NanoparticleAtomCounter.cli.atom_count",
                    "--input",
                    in_path,
                    "--output",
                    out_path,
                    "--mode",
                    mode,
                ],
                check=True,
                capture_output=True,
                text=True,
            )

        except subprocess.CalledProcessError as err:
            st.error(f"Atom-counter failed:\n{err.stderr}")
            os.remove(in_path)
            st.stop()

        # now, read and show results
        df_out = pd.read_csv(out_path)

        st.markdown("#### Results")
        st.download_button(
            "Download CSV",
            data=open(out_path, "rb").read(),
            file_name="atom_count_output.csv",
            mime="text/csv",
        )
        st.dataframe(df_out, use_container_width=True)

        # cleanup
        os.remove(in_path)
        os.remove(out_path)
