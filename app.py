import contextlib
import os
import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st
from nanoparticleatomcounting.atom_count import main as atom_counter


# ───────────────────────────  PAGE CONFIG  ───────────────────────────
st.set_page_config(page_title="Nanoparticle Atom Counter", page_icon="🧮")

# ---------- sidebar: resources ----------
with st.sidebar:
    st.header("Resources")

    SAMPLE_CSV = Path("sample_input.csv").read_bytes()
#    SAMPLE_CSV = (
#        "r (A),R (A),Theta,Element,Facet\n"
#        "10,,70.0,Ag,\n"
#        "100,,90,Pd,\n"
#        "400,,125,Cu,\n"
#    ).encode()

    st.download_button(
        "Sample input (.csv)",
        SAMPLE_CSV,
        file_name="sample_input.csv",
        mime="text/csv",
    )

    st.image("Acute.png", caption="θ < 90° (acute)", use_container_width=True)
    st.image("Obtuse.png", caption="θ > 90° (obtuse)", use_container_width=True)
    st.image("Nanoparticle_Legend.tif", caption="Definition of surface, interfacial, and perimeter atoms", use_container_width=True)
# ---------------------------------------------------------------------


# ───────────────────────────  MAIN LAYOUT  ────────────────────────────
st.title("Nanoparticle Atom Counter")

st.markdown(
    """
Upload a **.csv**, **.xls**, or **.xlsx** file containing the columns  
`r (A), R (A), Theta, Element, Interface Facet, Surface Facet`.

*Supply either **r** or **R** (if both are present, **r** is used).  
Interface Facet and Surface Facet are optional; leave blank if unknown.*

**Need a template or a visual guide?**  
A sample input file and explanatory diagrams are available in the sidebar.


**This note here is for me: when the paper has been written, add a link to the paper.md and the CITATION.cff**
""",
    unsafe_allow_html=True,
)

# counting-mode selector (now in main area)
mode = st.radio("Counting mode", ("volume", "area"), horizontal=True)


# ──────────  FILE UPLOAD  ──────────
file = st.file_uploader(
    "Drag-and-drop or browse for your input file",
    type=("csv", "xls", "xlsx"),
    accept_multiple_files=False,
)

if file is None:
    st.stop()      # wait for user input

"""
# ──────────  CALCULATION  ──────────
if st.button("⚙️ Run calculation"):
    with st.spinner("Processing …"):

        # temporary input file
        in_suffix = Path(file.name).suffix
        with tempfile.NamedTemporaryFile(delete=False, suffix=in_suffix) as tin:
            tin.write(file.getbuffer())
            tin.flush()
            in_path = tin.name

        # temporary output file
        with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as tout:
            out_path = tout.name

        # run CLI
        success = True
        with contextlib.suppress(Exception):
            atom_counter(in_path, out_path, mode=mode)
        if not Path(out_path).exists():
            success = False

        if not success:
            st.error("Calculation failed. Please check your input and try again.")
            os.remove(in_path)
            os.remove(out_path)
            st.stop()

        # read and show results
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
    """

# ──────────  CALCULATION  ──────────
if st.button("⚙️ Run calculation"):
    with st.spinner("Processing …"):

        # write the uploaded file to disk
        in_suffix = Path(file.name).suffix           # .csv / .xls / .xlsx
        fd_in, in_path  = tempfile.mkstemp(suffix=in_suffix)
        os.write(fd_in, file.getbuffer())
        os.close(fd_in)

        # pick a *name* for the output file, but don’t create it yet
        _, out_path = tempfile.mkstemp(suffix=".csv")
        os.remove(out_path)                          # ensure it doesn’t exist

        # ---------- run the CLI and catch any error ----------
        try:
            atom_counter(in_path, out_path, mode=mode)
        except Exception as e:
            st.exception(e)                          # show full traceback
            Path(in_path).unlink(missing_ok=True)
            st.stop()

        # ---------- sanity-check the output ----------
        if not Path(out_path).exists() or Path(out_path).stat().st_size == 0:
            st.error(
                "The calculation finished, but the output file is empty. "
                "Check your input values or look at the traceback above."
            )
            Path(in_path).unlink(missing_ok=True)
            Path(out_path).unlink(missing_ok=True)
            st.stop()

        # ---------- load & display ----------
        try:
            df_out = pd.read_csv(out_path)
        except pd.errors.EmptyDataError:
            st.error(f"{out_path} exists but contains no rows.")
            Path(in_path).unlink(missing_ok=True)
            Path(out_path).unlink(missing_ok=True)
            st.stop()

        st.markdown("#### Results")
        st.download_button(
            "Download CSV",
            data=Path(out_path).read_bytes(),
            file_name="atom_count_output.csv",
            mime="text/csv",
        )
        st.dataframe(df_out, use_container_width=True)

        # cleanup temp files
        Path(in_path).unlink(missing_ok=True)
        Path(out_path).unlink(missing_ok=True)

