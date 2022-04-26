import streamlit as st
from support.variant import variant
from support.ranking import ranking
from support.about import about

# Config the whole app
st.set_page_config(
    page_title="Diagnosis", page_icon="🧊", layout="wide", initial_sidebar_state="expanded",
)

st.write(
    "<style>div.row-widget.stRadio > div{flex-direction:row;justify-content: center;} </style>",
    unsafe_allow_html=True,
)
st.write(
    "<style>div.st-bf{flex-direction:column;} div.st-ag{font-weight:bold;padding-right:50px;}</style>",
    unsafe_allow_html=True,
)


def main():
    """Navigating streamlit app"""

    # st.sidebar.title("Tools")

    PAGES = {"Home": about, "Variant": variant, "Ranking": ranking}

    # Select pages
    # Use dropdown if you prefer
    selection = st.radio("", list(PAGES.keys()))

    page = PAGES[selection]

    with st.spinner(f"Loading Page {selection} ..."):
        page = page()


if __name__ == "__main__":
    main()
