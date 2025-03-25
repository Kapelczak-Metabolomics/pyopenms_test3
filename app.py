# app.py  
import streamlit as st  
import plotly.graph_objects as go  
from pyopenms import MSExperiment, MzMLFile  
import numpy as np  
import pandas as pd  
import io  
  
# Configure page  
st.set_page_config(page_title="PyOpenMS Mass Spectrometry Analyzer", layout="wide")  
st.title("Mass Spectrometry Data Analysis")  
  
# Helper Functions  
def load_mzml(uploaded_file):  
    experiment = MSExperiment()  
    mzml_file = MzMLFile()  
    try:  
        mzml_file.load(uploaded_file, experiment)  
        return experiment  
    except Exception as e:  
        st.error(f"Error loading mzML file: {str(e)}")  
        return None  
  
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
    mz_values = []  
    intensities = []  
      
    for line in io.StringIO(mgf_bytes.decode('utf-8')).readlines():  
        line = line.strip()  
        if line == "BEGIN IONS":  
            current_record = {}  
            mz_values = []  
            intensities = []  
        elif line == "END IONS":  
            if current_record is not None:  
                current_record['mz_array'] = np.array(mz_values)  
                current_record['intensity_array'] = np.array(intensities)  
                mgf_data.append(current_record)  
        elif "=" in line:  
            key, value = line.split("=", 1)  # Split on first equals sign only  
            key = key.strip().lower()  
            value = value.strip()  
            if key == "pepmass":  
                current_record[key] = float(value.split()[0])  
            elif key == "charge":  
                current_record[key] = value  
            else:  
                current_record[key] = value  
        elif line and current_record is not None:  
            try:  
                mz, intensity = map(float, line.split())  
                mz_values.append(mz)  
                intensities.append(intensity)  
            except ValueError:  
                pass  
    return mgf_data  
  
def match_fragments(ms2_spec, mgf_spec, tolerance=0.5):  
    ms2_indices = []  
    mgf_indices = []  
      
    for i, ms2_mz in enumerate(ms2_spec['mz']):  
        for j, mgf_mz in enumerate(mgf_spec['mz_array']):  
            if abs(ms2_mz - mgf_mz) <= tolerance:  
                ms2_indices.append(i)  
                mgf_indices.append(j)  
                break  
      
    return {  
        'ms2_indices': ms2_indices,  
        'mgf_indices': mgf_indices  
    }  
  
# Main app  
uploaded_file = st.file_uploader("Upload an mzML file", type=["mzml"])  
  
if uploaded_file is not None:  
    experiment = load_mzml(uploaded_file)  
      
    if experiment:  
        st.success(f"File loaded successfully: {uploaded_file.name}")  
          
        # Display TIC  
        st.header("Total Ion Chromatogram (TIC)")  
        times, intensities = extract_tic(experiment)  
          
        fig_tic = go.Figure()  
        fig_tic.add_trace(go.Scatter(  
            x=times,   
            y=intensities,  
            mode='lines',  
            line=dict(color="#2563EB", width=2)  
        ))  
        fig_tic.update_layout(  
            title={"text": "Total Ion Chromatogram", "pad": {"t":15}},  
            xaxis_title="Retention Time (s)",  
            yaxis_title="Intensity",  
            xaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
            yaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
            template="simple_white"  
        )  
        st.plotly_chart(fig_tic, use_container_width=True)  
          
        # XIC Extraction  
        st.header("Extracted Ion Chromatogram (XIC)")  
        col1, col2 = st.columns(2)  
        with col1:  
            mass_input = st.text_input("Enter m/z value(s) to extract (comma-separated)", "")  
        with col2:  
            tolerance = st.number_input("Mass tolerance (Da)", min_value=0.01, value=0.5)  
          
        if mass_input:  
            masses = [float(m.strip()) for m in mass_input.split(",")]  
              
            fig_xic = go.Figure()  
            colors = ["#2563EB", "#24EB84", "#B2EB24", "#EB3424", "#D324EB"]  
              
            for i, mass in enumerate(masses):  
                xic_times, xic_intensities = extract_xic(experiment, mass, tolerance)  
                fig_xic.add_trace(go.Scatter(  
                    x=xic_times,  
                    y=xic_intensities,  
                    mode='lines',  
                    name=f"m/z {mass}",  
                    line=dict(color=colors[i % len(colors)], width=2)  
                ))  
              
            fig_xic.update_layout(  
                title={"text": "Extracted Ion Chromatogram", "pad": {"t":15}},  
                xaxis_title="Retention Time (s)",  
                yaxis_title="Intensity",  
                xaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
                yaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
                template="simple_white"  
            )  
            st.plotly_chart(fig_xic, use_container_width=True)  
          
        # MS2 Data  
        st.header("MS2 Fragmentation Data")  
        ms2_data = extract_ms2_data(experiment)  
          
        if ms2_data:  
            ms2_idx = st.selectbox("Select an MS2 spectrum", list(range(len(ms2_data))))  
            selected_ms2 = ms2_data[ms2_idx]  
            st.write(f"Precursor m/z: {selected_ms2['precursor']}")  
              
            fig_ms2 = go.Figure()  
            fig_ms2.add_trace(go.Bar(  
                x=selected_ms2['mz'],  
                y=selected_ms2['intensity'],  
                marker_color="#B2EB24"  
            ))  
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
        uploaded_mgf = st.file_uploader("Upload an MGF file", type=["mgf"])  
          
        if uploaded_mgf is not None:  
            try:  
                mgf_records = parse_mgf(uploaded_mgf.getvalue())  
                st.success(f"MGF file parsed successfully: {len(mgf_records)} records found")  
                  
                record_idx = st.selectbox("Select an MGF record", list(range(len(mgf_records))))  
                record = mgf_records[record_idx]  
                  
                st.write(f"Peptide Mass: {record.get('pepmass', 'N/A')}")  
                st.write(f"Charge: {record.get('charge', 'N/A')}")  
                  
                fig_mgf = go.Figure()  
                fig_mgf.add_trace(go.Bar(  
                    x=record['mz_array'],  
                    y=record['intensity_array'],  
                    marker_color="#EB3424"  
                ))  
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
                    match_tol = st.number_input("Matching tolerance (Da)", min_value=0.01, value=0.5, key="match_tol")  
                    ms2_match_idx = st.selectbox("Select MS2 spectrum for matching", list(range(len(ms2_data))), key="ms2_match")  
                      
                    ms2_spec = ms2_data[ms2_match_idx]  
                    match_results = match_fragments(ms2_spec, record, match_tol)  
                      
                    if match_results['ms2_indices']:  
                        st.success(f"Found {len(match_results['ms2_indices'])} matching fragments")  
                          
                        # Create visualization of matches  
                        fig_match = go.Figure()  
                          
                        # MS2 spectrum (positive y-axis)  
                        fig_match.add_trace(go.Bar(  
                            x=ms2_spec['mz'],  
                            y=ms2_spec['intensity'],  
                            name="MS2 Spectrum",  
                            marker_color="#B2EB24"  
                        ))  
                          
                        # MGF spectrum (negative y-axis for comparison)  
                        fig_match.add_trace(go.Bar(  
                            x=record['mz_array'],  
                            y=-record['intensity_array'],  
                            name="MGF Spectrum",  
                            marker_color="#EB3424"  
                        ))  
                          
                        # Add connecting lines for matches  
                        for ms2_idx, mgf_idx in zip(match_results['ms2_indices'], match_results['mgf_indices']):  
                            ms2_mz = ms2_spec['mz'][ms2_idx]  
                            ms2_intensity = ms2_spec['intensity'][ms2_idx]  
                            mgf_mz = record['mz_array'][mgf_idx]  
                            mgf_intensity = -record['intensity_array'][mgf_idx]  
                              
                            fig_match.add_shape(  
                                type="line",  
                                x0=ms2_mz, y0=ms2_intensity,  
                                x1=mgf_mz, y1=mgf_intensity,  
                                line=dict(color="#2563EB", width=1, dash="dot")  
                            )  
                          
                        fig_match.update_layout(  
                            title={"text": "Fragment Matching Visualization", "pad": {"t":15}},  
                            xaxis_title="m/z",  
                            yaxis_title="Intensity",  
                            xaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
                            yaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
                            template="simple_white"  
                        )  
                        st.plotly_chart(fig_match, use_container_width=True)  
                          
                        # Show matching table  
                        match_data = {  
                            "MS2 m/z": [ms2_spec['mz'][i] for i in match_results['ms2_indices']],  
                            "MGF m/z": [record['mz_array'][j] for j in match_results['mgf_indices']],  
                            "Difference": [ms2_spec['mz'][i] - record['mz_array'][j] for i, j in zip(match_results['ms2_indices'], match_results['mgf_indices'])]  
                        }  
                        match_df = pd.DataFrame(match_data)  
                        st.dataframe(match_df)  
                    else:  
                        st.warning("No matching fragments found with the current tolerance.")  
            except Exception as e:  
                st.error(f"Error processing MGF file: {str(e)}")  
else:  
    st.info("Please upload an mzML file to begin analysis.")  
      
    # Show features  
    st.subheader("Features")  
    st.markdown("""  
    - **Total Ion Chromatogram (TIC)**: View the total ion current across retention time  
    - **Extracted Ion Chromatogram (XIC)**: Extract and visualize specific m/z values  
    - **MS2 Fragmentation Analysis**: View MS2 spectra and analyze fragmentation patterns  
    - **MGF Matching**: Match experimental MS2 data with theoretical fragments from MGF files  
    """)  
