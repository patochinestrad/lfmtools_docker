import streamlit as st
import json


with open("config.json", "r") as configfile:
    configurations = json.load(configfile)

for i in configurations:
    if i not in st.session_state:
        st.session_state[i] = configurations[i]


st.title("LFM DrugDesign App")
st.subheader(
    "To change your app configuration modify the config.json file. Current configuration: "
)


for i in configurations:
    c = st.container()
    with c:
        st.header("%s path" % i)
        st.code(configurations[i])
