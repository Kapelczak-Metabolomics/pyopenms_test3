# app.py  
import streamlit as st  
import plotly.graph_objects as go  
from pyopenms import MSExperiment, MzMLFile, MSSpectrum  
import numpy as np  
import pandas as pd  
  
# Configure page  
st.set_page_config(page_title="PyOpenMS Mass Spectrometry Analyzer", layout="wide")  
st.title("Mass Spectrometry Data Analysis")  
  
# Helper Functions  
def load_mzml(uploaded_file):  
    experiment = MSExperiment()  
    mzml_file = MzMLFile()  
    try:  
        mzml_file.load(uploaded_file, experiment)  
    except Exception as e:  
        st.error(f"Error loading mzML file: {str(e)}")  
    return experiment  
  
def extract_tic(experiment):  
    times, intensities = [], []  
    for spectrum in experiment:  
        if spectrum.getMSLevel() == 1:  
            time = spectrum.getRT()  
            peaks = spectrum.get_peaks()  
            intensity = sum(peaks[1]) if peaks[1].size > 0 else 0.0  
            times.append(time)  
            intensities.append(intensity)  
    return np.array(times), np.array(intensities)  
  
def extract_xic(experiment, mass, tol=0.5):  
    times, intensities = [], []  
    for spectrum in experiment:  
        if spectrum.getMSLevel() == 1:  
            time = spectrum.getRT()  
            mz_array, intensity_array = spectrum.get_peaks()  
            mask = (mz_array >= mass - tol) & (mz_array <= mass + tol)  
            intensity = np.sum(intensity_array[mask]) if np.any(mask) else 0.0  
            times.append(time)  
            intensities.append(intensity)  
    return np.array(times), np.array(intensities)  
  
def extract_ms2_data(experiment):  
    ms2_data = []  
    for i, spectrum in enumerate(experiment):  
        if spectrum.getMSLevel() == 2:  
            precursors = spectrum.getPrecursors()  
            if precursors:  
                precursor_mz = precursors[0].getMZ()  
                mz_array, intensity_array = spectrum.get_peaks()  
                ms2_data.append({  
                    'index': i,  
                    'rt': spectrum.getRT(),  
                    'precursor': precursor_mz,  
                    'mz': mz_array,  
                    'intensity': intensity_array  
                })  
    return ms2_data  
  
def parse_mgf(mgf_bytes):  
    mgf_data = []  
    current_record = None  
    mz_values, intensities = [], []  
      
    for line in mgf_bytes.decode('utf-8').splitlines():  
        line = line.strip()  
        if line == "BEGIN IONS":  
            current_record = {}  
            mz_values, intensities = [], []  
        elif line == "END IONS":  
            if current_record is not None:  
                current_record['mz_array'] = np.array(mz_values)  
                current_record['intensity_array'] = np.array(intensities)  
                mgf_data.append(current_record)  
        elif "=" in line:  
            try:  
                key, value = line.split("=", 1)  
                key = key.strip().lower()  
                value = value.strip()  
                if key == "pepmass":  
                    current_record[key] = float(value.split()[0])  
                elif key == "charge":  
                    current_record[key] = int(value.replace('+', '').replace('-', ''))  
                else:  
                    current_record[key] = value  
            except Exception:  
                pass  
        elif current_record is not None and line and not line.startswith("#"):  
            try:  
                parts = line.split()  
                if len(parts) >= 2:  
                    mz, intensity = float(parts[0]), float(parts[1])  
                    mz_values.append(mz)  
                    intensities.append(intensity)  
            except Exception:  
                pass  
      
    return mgf_data  
  
def match_fragments(ms2_spec, mgf_spec, tolerance=0.5):  
    matches = []  
    ms2_indices = []  
    mgf_indices = []  
      
    for i, mz1 in enumerate(ms2_spec['mz']):  
        for j, mz2 in enumerate(mgf_spec['mz_array']):  
            if abs(mz1 - mz2) <= tolerance:  
                matches.append((mz1, mz2))  
                ms2_indices.append(i)  
                mgf_indices.append(j)  
                break  
      
    return {  
        'matches': matches,  
        'ms2_indices': ms2_indices,  
        'mgf_indices': mgf_indices  
    }  
  
# Main app  
uploaded_file = st.file_uploader("Upload an mzML file", type=["mzml"])  
  
if uploaded_file is not None:  
    experiment = load_mzml(uploaded_file)  
    st.success(f"File loaded successfully: {uploaded_file.name}")  
      
    # TIC Display  
    st.header("Total Ion Chromatogram (TIC)")  
    times, intensities = extract_tic(experiment)  
      
    fig_tic = go.Figure()  
    fig_tic.add_trace(go.Scatter(x=times, y=intensities, mode='lines', line=dict(color="#2563EB", width=2)))  
    fig_tic.update_layout(  
        title={"text": "Total Ion Chromatogram", "pad": {"t":15}},  
        xaxis_title="Retention Time (min)",  
        yaxis_title="Intensity",  
        xaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
        yaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
        template="simple_white"  
    )  
    st.plotly_chart(fig_tic, use_container_width=True)  
      
    # XIC Extraction  
    st.header("Extracted Ion Chromatogram (XIC)")  
    col1, col2 = st.columns([3, 1])  
    with col1:  
        masses_input = st.text_input("Enter m/z values (comma-separated)", "")  
    with col2:  
        tolerance = st.number_input("m/z Tolerance (Da)", min_value=0.01, max_value=10.0, value=0.5, step=0.1)  
      
    if masses_input:  
        try:  
            masses = [float(m.strip()) for m in masses_input.split(",")]  
            fig_xic = go.Figure()  
              
            colors = ["#2563EB", "#24EB84", "#B2EB24", "#EB3424", "#D324EB"]  
              
            for i, mass in enumerate(masses):  
                xic_times, xic_intensities = extract_xic(experiment, mass, tolerance)  
                color = colors[i % len(colors)]  
                fig_xic.add_trace(go.Scatter(  
                    x=xic_times,   
                    y=xic_intensities,   
                    mode='lines',   
                    name=f"m/z {mass}",  
                    line=dict(color=color, width=2)  
                ))  
              
            fig_xic.update_layout(  
                title={"text": "Extracted Ion Chromatogram", "pad": {"t":15}},  
                xaxis_title="Retention Time (min)",  
                yaxis_title="Intensity",  
                xaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
                yaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
                template="simple_white"  
            )  
            st.plotly_chart(fig_xic, use_container_width=True)  
        except Exception as e:  
            st.error(f"Error extracting XIC: {str(e)}")  
      
    # MS2 Data  
    st.header("MS2 Fragmentation Data")  
    ms2_data = extract_ms2_data(experiment)  
      
    if ms2_data:  
        ms2_idx = st.selectbox("Select an MS2 spectrum", list(range(len(ms2_data))))  
        selected_ms2 = ms2_data[ms2_idx]  
        st.write(f"Precursor m/z: {selected_ms2['precursor']}")  
        st.write(f"Retention Time: {selected_ms2['rt']} min")  
          
        fig_ms2 = go.Figure()  
        fig_ms2.add_trace(go.Bar(x=selected_ms2['mz'], y=selected_ms2['intensity'], marker_color="#B2EB24"))  
        fig_ms2.update_layout(  
            title={"text": "MS2 Fragmentation Spectrum", "pad": {"t":15}},  
            xaxis_title="m/z",  
            yaxis_title="Intensity",  
            xaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
            yaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
            template="simple_white"  
        )  
        st.plotly_chart(fig_ms2, use_container_width=True)  
    else:  
        st.info("No MS2 data found in the uploaded mzML file.")  
      
    # MGF Matching  
    st.header("MGF Fragmentation Matching")  
    uploaded_mgf = st.file_uploader("Upload an MGF file for fragmentation matching", type=["mgf"])  
      
    if uploaded_mgf is not None:  
        try:  
            mgf_records = parse_mgf(uploaded_mgf.getvalue())  
            st.success(f"MGF file parsed successfully: {len(mgf_records)} records found")  
              
            if mgf_records:  
                record_idx = st.selectbox("Select an MGF record", list(range(len(mgf_records))))  
                record = mgf_records[record_idx]  
                  
                st.write(f"Peptide Mass: {record.get('pepmass', 'N/A')}")  
                st.write(f"Charge: {record.get('charge', 'N/A')}")  
                  
                fig_mgf = go.Figure()  
                fig_mgf.add_trace(go.Bar(x=record['mz_array'], y=record['intensity_array'], marker_color="#EB3424"))  
                fig_mgf.update_layout(  
                    title={"text": "MGF Fragmentation Pattern", "pad": {"t":15}},  
                    xaxis_title="m/z",  
                    yaxis_title="Intensity",  
                    xaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
                    yaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
                    template="simple_white"  
                )  
                st.plotly_chart(fig_mgf, use_container_width=True)  
                  
                # Matching  
                if ms2_data:  
                    st.subheader("Fragment Matching")  
                    match_ms2_idx = st.selectbox("Select MS2 spectrum for matching", list(range(len(ms2_data))), key="match_ms2")  
                    ms2_spec = ms2_data[match_ms2_idx]  
                      
                    match_tolerance = st.number_input("Matching Tolerance (Da)", min_value=0.01, max_value=2.0, value=0.5, step=0.1)  
                    match_results = match_fragments(ms2_spec, record, match_tolerance)  
                      
                    if match_results['matches']:  
                        st.success(f"Found {len(match_results['matches'])} matching fragments")  
                          
                        # Create a combined plot  
                        fig_match = go.Figure()  
                          
                        # MS2 data  
                        fig_match.add_trace(go.Bar(  
                            x=ms2_spec['mz'],   
                            y=ms2_spec['intensity'],  
                            name="MS2 Spectrum",  
                            marker_color="#B2EB24",  
                            opacity=0.7  
                        ))  
                          
                        # MGF data  
                        fig_match.add_trace(go.Bar(  
                            x=record['mz_array'],   
                            y=-record['intensity_array'],  # Negative to show below x-axis  
                            name="MGF Spectrum",  
                            marker_color="#EB3424",  
                            opacity=0.7  
                        ))  
                          
                        # Highlight matches  
                        for ms2_idx, mgf_idx in zip(match_results['ms2_indices'], match_results['mgf_indices']):  
                            ms2_mz = ms2_spec['mz'][ms2_idx]  
                            ms2_intensity = ms2_spec['intensity'][ms2_idx]  
                            mgf_mz = record['mz_array'][mgf_idx]  
                            mgf_intensity = record['intensity_array'][mgf_idx]  
                              
                            # Add vertical line connecting matches  
                            fig_match.add_shape(  
                                type="line",  
                                x0=ms2_mz, y0=ms2_intensity,  
                                x1=mgf_mz, y1=-mgf_intensity,  
                                line=dict(color="#2563EB", width=1, dash="dot")  
                            )  
                          
                        fig_match.update_layout(  
                            title={"text": "Fragment Matching Visualization", "pad": {"t":15}},  
                            xaxis_title="m/z",  
                            yaxis_title="Intensity",  
                            xaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
                            yaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
                            template="simple_white",  
                            barmode='overlay'  
                        )  
                        st.plotly_chart(fig_match, use_container_width=True)  
                ▍
