# app.py  
import streamlit as st  
import plotly.graph_objects as go  
from pyopenms import MSExperiment, MzMLFile, MSSpectrum  
import numpy as np  
import io  
import re  
  
st.set_page_config(page_title="PyOpenMS Mass Spectrometry Analyzer", layout="wide")  
st.title("Mass Spectrometry Data Analysis")  
st.write("Upload an mzML file to analyze mass spectrometry data and perform XIC extraction")  
  
# --- Helper Functions ---  
def load_mzml(file_bytes):  
    """  
    Load an mzML file from uploaded bytes and return an MSExperiment.  
    """  
    mzml_data = MSExperiment()  
    mzml_file = MzMLFile()  
    mzml_file.loadFromMemory(file_bytes.read(), mzml_data)  
    return mzml_data  
  
def extract_tic(experiment):  
    """  
    Extract the Total Ion Chromatogram (TIC) from an MSExperiment.  
    Returns times and intensities.  
    """  
    times = []  
    intensities = []  
    # Iterate over individual spectra  
    for spectrum in experiment:  
        if spectrum.getMSLevel() == 1:  
            time = spectrum.getRT()  # retention time  
            intens = sum(spectrum.get_peaks()[1]) if len(spectrum.get_peaks()[1]) > 0 else 0.0  
            times.append(time)  
            intensities.append(intens)  
    return np.array(times), np.array(intensities)  
  
def extract_xic(experiment, mass, tol=0.5):  
    """  
    Extract an extracted ion chromatogram (XIC) for a given mass.  
    The tolerance is in Da.  
    Returns times and intensities.  
    """  
    times = []  
    intensities = []  
    for spectrum in experiment:  
        # Only consider MS1 spectra for XIC extraction.  
        if spectrum.getMSLevel() == 1:  
            time = spectrum.getRT()  
            mz_array, intensity_array = spectrum.get_peaks()  
            mz_array = np.array(mz_array)  
            intensity_array = np.array(intensity_array)  
            # Select intensities within the tolerance window  
            selection = (mz_array >= mass - tol) & (mz_array <= mass + tol)  
            intens = np.sum(intensity_array[selection]) if np.any(selection) else 0.0  
            times.append(time)  
            intensities.append(intens)  
    return np.array(times), np.array(intensities)  
  
def get_ms2_data(experiment):  
    """  
    Extract MS2 data from an MSExperiment.  
    Returns a list of dictionaries with MS2 data.  
    """  
    ms2_data = []  
    for spectrum in experiment:  
        if spectrum.getMSLevel() == 2:  
            precursors = spectrum.getPrecursors()  
            if precursors:  
                precursor_mz = precursors[0].getMZ()  
                mz_array, intensity_array = spectrum.get_peaks()  
                ms2_data.append({  
                    'rt': spectrum.getRT(),  
                    'precursor': precursor_mz,  
                    'mz': mz_array,  
                    'intensity': intensity_array  
                })  
    return ms2_data  
  
def parse_mgf(mgf_bytes):  
    """  
    Parse an MGF file and return a list of dictionaries with MGF data.  
    Each dictionary represents a spectrum with its metadata and peak list.  
    """  
    mgf_data = []  
    current_record = None  
    mz_values = []  
    intensities = []  
      
    # Decode bytes to string and split by lines  
    mgf_text = mgf_bytes.decode('utf-8')  
    lines = mgf_text.split('\n')  
      
    for line in lines:  
        line = line.strip()  
        if not line:  
            continue  
              
        if line == "BEGIN IONS":  
            current_record = {}  
            mz_values = []  
            intensities = []  
        elif line == "END IONS":  
            if current_record is not None:  
                current_record['mz_array'] = np.array(mz_values)  
                current_record['intensity_array'] = np.array(intensities)  
                mgf_data.append(current_record)  
                current_record = None  
        elif "=" in line:  
            # Fix: Split only on the first equals sign  
            parts = line.split("=", 1)  
            if len(parts) == 2:  
                key, value = parts  
                key = key.strip().lower()  
                value = value.strip()  
                  
                # Convert appropriate types  
                if key == "pepmass":  
                    # Handle pepmass which might have intensity after space  
                    pepmass_parts = value.split()  
                    value = float(pepmass_parts[0])  
                elif key == "charge":  
                    # Remove + or - and convert to int  
                    value = int(re.sub(r'[+-]', '', value))  
                  
                current_record[key] = value  
        elif current_record is not None and " " in line:  
            # This is a peak line (m/z intensity)  
            try:  
                parts = line.split()  
                if len(parts) >= 2:  
                    mz = float(parts[0])  
                    intensity = float(parts[1])  
                    mz_values.append(mz)  
                    intensities.append(intensity)  
            except ValueError:  
                # Skip lines that can't be parsed as peaks  
                pass  
      
    return mgf_data  
  
def match_fragments(ms2_spectrum, mgf_spectrum, tolerance=0.5):  
    """  
    Match fragments between MS2 spectrum and MGF spectrum.  
    Returns matched peaks and their intensities.  
    """  
    ms2_mz = ms2_spectrum['mz']  
    mgf_mz = mgf_spectrum['mz_array']  
      
    matched_indices = []  
    matched_mgf_indices = []  
      
    for i, mz in enumerate(ms2_mz):  
        # Find closest match in MGF spectrum  
        diffs = np.abs(mgf_mz - mz)  
        min_idx = np.argmin(diffs)  
          
        if diffs[min_idx] <= tolerance:  
            matched_indices.append(i)  
            matched_mgf_indices.append(min_idx)  
      
    return {  
        'ms2_indices': matched_indices,  
        'mgf_indices': matched_mgf_indices,  
        'match_count': len(matched_indices),  
        'total_ms2_peaks': len(ms2_mz),  
        'total_mgf_peaks': len(mgf_mz)  
    }  
  
# Main app  
st.sidebar.header("Settings")  
  
# File uploader for mzML file  
uploaded_file = st.file_uploader("Choose an mzML file", type=["mzml"])  
  
if uploaded_file is not None:  
    # Load the mzML file  
    with st.spinner("Loading mzML file..."):  
        try:  
            experiment = load_mzml(uploaded_file)  
            st.success(f"Successfully loaded mzML file with {experiment.size()} spectra")  
              
            # Extract and display TIC  
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
                title={"text": "Total Ion Chromatogram", "font": {"size": 20, "family": "Inter"}, "pad": {"t": 15}},  
                xaxis_title={"text": "Retention Time (seconds)", "font": {"size": 16, "family": "Inter"}, "standoff": 10},  
                yaxis_title={"text": "Intensity", "font": {"size": 16, "family": "Inter"}, "standoff": 10},  
                xaxis=dict(  
                    showgrid=False,   
                    gridcolor="#F3F4F6",   
                    showline=True,   
                    linecolor="#171717",   
                    linewidth=1,  
                    ticks="outside",  
                    tickfont=dict(family="Inter", size=14, color="#171717")  
                ),  
                yaxis=dict(  
                    showgrid=True,   
                    gridcolor="#F3F4F6",   
                    showline=True,   
                    linecolor="#171717",   
                    linewidth=1,  
                    ticks="outside",  
                    tickfont=dict(family="Inter", size=14, color="#171717")  
                ),  
                plot_bgcolor='white',  
                margin=dict(l=10, r=10, t=50, b=10),  
            )  
            # Hide top and right spines  
            fig_tic.update_xaxes(showspikes=True, mirror=False, showline=True, linewidth=1, linecolor="#171717")  
            fig_tic.update_yaxes(showspikes=True, mirror=False, showline=True, linewidth=1, linecolor="#171717")  
              
            st.plotly_chart(fig_tic, use_container_width=True)  
              
            # XIC Extraction  
            st.header("Extracted Ion Chromatogram (XIC)")  
            col1, col2 = st.columns([3, 1])  
              
            with col1:  
                mass_input = st.text_input("Enter m/z values (comma-separated)", "")  
              
            with col2:  
                tolerance = st.number_input("Mass Tolerance (Da)", min_value=0.01, max_value=10.0, value=0.5, step=0.1)  
              
            if mass_input:  
                masses = [float(m.strip()) for m in mass_input.split(",") if m.strip()]  
                  
                if masses:  
                    fig_xic = go.Figure()  
                      
                    # Color palette for multiple XICs  
                    colors = ["#2563EB", "#24EB84", "#B2EB24", "#EB3424", "#D324EB"]  
                      
                    for i, mass in enumerate(masses):  
                        xic_times, xic_intensities = extract_xic(experiment, mass, tolerance)  
                        color_idx = i % len(colors)  
                          
                        fig_xic.add_trace(go.Scatter(  
                            x=xic_times,   
                            y=xic_intensities,   
                            mode='lines',   
                            name=f"m/z {mass:.4f}",  
                            line=dict(color=colors[color_idx], width=2)  
                        ))  
                      
                    fig_xic.update_layout(  
                        title={"text": "Extracted Ion Chromatogram", "font": {"size": 20, "family": "Inter"}, "pad": {"t": 15}},  
                        xaxis_title={"text": "Retention Time (seconds)", "font": {"size": 16, "family": "Inter"}, "standoff": 10},  
                        yaxis_title={"text": "Intensity", "font": {"size": 16, "family": "Inter"}, "standoff": 10},  
                        xaxis=dict(  
                            showgrid=False,   
                            gridcolor="#F3F4F6",   
                            showline=True,   
                            linecolor="#171717",   
                            linewidth=1,  
                            ticks="outside",  
                            tickfont=dict(family="Inter", size=14, color="#171717")  
                        ),  
                        yaxis=dict(  
                            showgrid=True,   
                            gridcolor="#F3F4F6",   
                            showline=True,   
                            linecolor="#171717",   
                            linewidth=1,  
                            ticks="outside",  
                            tickfont=dict(family="Inter", size=14, color="#171717")  
                        ),  
                        plot_bgcolor='white',  
                        margin=dict(l=10, r=10, t=50, b=10),  
                        legend=dict(  
                            font=dict(family="Inter", size=14, color="#171717"),  
                            orientation="h",  
                            yanchor="bottom",  
                            y=1.02,  
                            xanchor="right",  
                            x=1  
                        )  
                    )  
                    # Hide top and right spines  
                    fig_xic.update_xaxes(showspikes=True, mirror=False, showline=True, linewidth=1, linecolor="#171717")  
                    fig_xic.update_yaxes(showspikes=True, mirror=False, showline=True, linewidth=1, linecolor="#171717")  
                      
                    st.plotly_chart(fig_xic, use_container_width=True)  
              
            # MS2 Data  
            st.header("MS2 Fragmentation Data")  
            ms2_data = get_ms2_data(experiment)  
              
            if ms2_data:  
                st.write(f"Found {len(ms2_data)} MS2 spectra")  
                  
                # Create a selectbox with precursor m/z values  
                precursor_options = [f"RT: {spec['rt']:.2f} - Precursor m/z: {spec['precursor']:.4f}" for spec in ms2_data]  
                ms2_idx = st.selectbox("Select MS2 spectrum", range(len(precursor_options)), format_func=lambda i: precursor_options[i])  
                  
                selected_ms2 = ms2_data[ms2_idx]  
                  
                # Display MS2 spectrum  
                fig_ms2 = go.Figure()  
                fig_ms2.add_trace(go.Bar(  
                    x=selected_ms2['mz'],   
                    y=selected_ms2['intensity'],   
                    marker_color="#B2EB24"  
                ))  
                  
                fig_ms2.update_layout(  
                    title={"text": "MS2 Fragmentation Spectrum", "font": {"size": 20, "family": "Inter"}, "pad": {"t": 15}},  
                    xaxis_title={"text": "m/z", "font": {"size": 16, "family": "Inter"}, "standoff": 10},  
                    yaxis_title={"text": "Intensity", "font": {"size": 16, "family": "Inter"}, "standoff": 10},  
                    xaxis=dict(  
                        showgrid=False,   
                        gridcolor="#F3F4F6",   
                        showline=True,   
                        linecolor="#171717",   
                        linewidth=1,  
                        ticks="outside",  
                        tickfont=dict(family="Inter", size=14, color="#171717")  
                    ),  
                    yaxis=dict(  
                        showgrid=True,   
                        gridcolor="#F3F4F6",   
                        showline=True,   
                        linecolor="#171717",   
                        linewidth=1,  
                        ticks="outside",  
                        tickfont=dict(family="Inter", size=14, color="#171717")  
                    ),  
                    plot_bgcolor='white',  
                    margin=dict(l=10, r=10, t=50, b=10),  
                )  
                # Hide top and right spines  
                fig_ms2.update_xaxes(showspikes=True, mirror=False, showline=True, linewidth=1, linecolor="#171717")  
                fig_ms2.update_yaxes(showspikes=True, mirror=False, showline=True, linewidth=1, linecolor="#171717")  
                  
                st.plotly_chart(fig_ms2, use_container_width=True)  
                  
                # Display peak table  
                with st.expander("Show MS2 Peak Table"):  
                    peak_data = {  
                        "m/z": selected_ms2['mz'],  
                        "Intensity": selected_ms2['intensity']  
                    }  
                    peak_df = pd.DataFrame(peak_data)  
                    st.dataframe(peak_df)  
            else:  
                st.info("No MS2 data found in the uploaded mzML file.")  
  
            # File uploader for MGF file  
            st.header("MGF Fragmentation Matching")  
            uploaded_mgf = st.file_uploader("Choose an MGF file for fragmentation matching", type=["mgf"])  
              
            if uploaded_mgf is not None:  
                with st.spinner("Parsing MGF file for fragmentation matching..."):  
                    try:  
                        mgf_records = parse_mgf(uploaded_mgf.getvalue())  
                        st.success(f"MGF file parsed successfully! Found {len(mgf_records)} spectra.")  
                          
                        if mgf_records and ms2_data:  
                            st.subheader("MGF Fragmentation Data")  
                              
                            # Create a selectbox for MGF records  
                            mgf_options = [f"PepMass: {rec.get('pepmass', 'N/A')}" for rec in mgf_records]  
                            mgf_idx = st.selectbox("Select an MGF record", range(len(mgf_options)), format_func=lambda i: mgf_options[i])  
                              
                            record = mgf_records[mgf_idx]  
                              
                            # Display metadata  
                            col1, col2 = st.columns(2)  
                            with col1:  
                                st.write(f"Peptide Mass: {record.get('pepmass', 'N/A')}")  
                            with col2:  
                                st.write(f"Charge: {record.get('charge', 'N/A')}")  
                              
                            # Display fragmentation pattern  
                            fig_mgf = go.Figure()  
                            fig_mgf.add_trace(go.Bar(  
                                x=record['mz_array'],   
                                y=record['intensity_array'],   
                                marker_color="#EB3424"  
                            ))  
                              
                            fig_mgf.update_layout(  
                                title={"text": "MGF Fragmentation Pattern", "font": {"size": 20, "family": "Inter"}, "pad": {"t": 15}},  
                                xaxis_title={"text": "m/z", "font": {"size": 16, "family": "Inter"}, "standoff": 10},  
                                yaxis_title={"text": "Intensity", "font": {"size": 16, "family": "Inter"}, "standoff": 10},  
                                xaxis=dict(  
                                    showgrid=False,   
                                    gridcolor="#F3F4F6",   
                                    showline=True,   
                                    linecolor="#171717",   
                                    linewidth=1,  
                                    ticks="outside",  
                                    tickfont=dict(family="Inter", size=14, color="#171717")  
                                ),  
                                yaxis=dict(  
                                    showgrid=True,   
                                    gridcolor="#F3F4F6",   
                                    showline=True,   
                                    linecolor="#171717",   
                                    linewidth=1,  
                                    ticks="outside",  
                                    tickfont=dict(family="Inter", size=14, color="#171717")  
                                ),  
                                plot_bgcolor='white',  
                                margin=dict(l=10, r=10, t=50, b=10),  
                            )  
                            # Hide top and right spines  
                            fig_mgf.update_xaxes(showspikes=True, mirror=False, showline=True, linewidth=1, linecolor="#171717")  
                            fig_mgf.update_yaxes(showspikes=True, mirror=False, showline=True, linewidth=1, linecolor="#171717")  
                              
                            st.plotly_chart(fig_mgf, use_container_width=True)  
                              
                            # Fragmentation Matching  
                            st.subheader("Fragmentation Matching")  
                              
                            # Allow user to select MS2 spectrum for matching  
                            match_ms2_idx = st.selectbox(  
                                "Select MS2 spectrum for matching",   
                                range(len(precursor_options)),   
                                format_func=lambda i: precursor_options[i],  
                                key="match_ms2"  
                            )  
                              
                            # Allow user to select MGF spectrum for matching  
                            match_mgf_idx = st.selectbox(  
                                "Select MGF spectrum for matching",   
                                range(len(mgf_options)),   
                                format_func=lambda i: mgf_options[i],  
                                key="match_mgf"  
                            )  
                              
                            # Set tolerance for matching  
                            match_tolerance = st.slider(  
                                "Matching Tolerance (Da)",   
                                min_value=0.01,   
                                max_value=2.0,   
                                value=0.5,   
                                step=0.01  
                            )  
                              
                            # Perform matching  
                            if st.button("Match Fragments"):  
                                ms2_spec = ms2_data[match_ms2_idx]  
                                mgf_spec = mgf_records[match_mgf_idx]  
                                  
                                match_results = match_fragments(ms2_spec, mgf_spec, match_tolerance)  
                                  
                                # Display match statistics  
                                st.write(f"Matched {match_results['match_count']} out of {match_results['total_ms2_peaks']} MS2 peaks")  
                                st.write(f"Match percentage: {match_results['match_count'] / match_results['total_ms2_peaks'] * 100:.2f}%")  
                                  
                                # Create a combined plot showing both spectra and matches  
                                fig_match = go.Figure()  
                                  
                                # Add MS2 spectrum  
                                fig_match.add_trace(go.Bar(  
                                    x=ms2_spec['mz'],  
                                    y=ms2_spec['intensity'],  
                                    name="MS2 Spectrum",  
                                    marker_color="#B2EB24",  
                                    opacity=0.7  
                                ))  
                                  
                                # Add MGF spectrum (inverted for clarity)  
                                mgf_intensities_inverted = -1 * mgf_spec['intensity_array']  
                                fig_match.add_trace(go.Bar(  
                                    x=mgf_spec['mz_array'],  
                                    y=mgf_intensities_inverted,  
                                    name="MGF Spectrum",  
                                    marker_color="#EB3424",  
                                    opacity=0.7  
                                ))  
                                  
                                # Highlight matched peaks in MS2  
                                if match_results['ms2_indices']:  
                                    matched_mz = [ms2_spec['mz'][i] for i in match_results['ms2_indices']]  
                                    matched_intensity = [ms2_spec['intensity'][i] for i in match_results['ms2_indices']]  
                                      
                                    fig_match.add_trace(go.Scatter(  
                                        x=matched_mz,  
                                        y=matched_intensity,  
                                        mode="markers",  
                                        name="Matched MS2 Peaks",  
                                        marker=dict(  
                                            color="#D324EB",  
                                            size=10,  
                                            symbol="circle"  
                                        )  
                                    ))  
                                  
                                # Highlight matched peaks in MGF  
                                if match_results['mgf_indices']:  
                                    matched_mgf_mz = [mgf_spec['mz_array'][i] for i in match_results['mgf_indices']]  
                                    matched_mgf_intensity = [-1 * mgf_spec['intensity_array'][i] for i in match_results['mgf_indices']]  
                                      
                                    fig_match.add_trace(go.Scatter(  
                                        x=matched_mgf_mz,  
                                        y=matched_mgf_intensity,  
                                        mode="markers",  
                                        name="Matched MGF Peaks",  
                                        marker=dict(  
                                            color="#D324EB",  
                                            size=10,  
                                            symbol="circle"  
                                        )  
                                    ))  
                                  
                                fig_match.update_layout(  
                                    title={"text": "Fragment Matching Visualization", "font": {"size": 20, "family": "Inter"}, "pad": {"t": 15}},  
                                    xaxis_title={"text": "m/z", "font": {"size": 16, "family": "Inter"}, "standoff": 10},  
                                    yaxis_title={"text": "Intensity (MGF inverted)", "font": {"size": 16, "family": "Inter"}, "standoff": 10},  
                                    xaxis=dict(  
                                        showgrid=False,   
                                        gridcolor="#F3F4F6",   
                                        showline=True,   
                                        linecolor="#171717",   
                                        linewidth=1,  
                                        ticks="outside",  
                                        tickfont=dict(family="Inter", size=14, color="#171717")  
                                    ),  
                                    yaxis=dict(  
                                        showgrid=True,   
                                        gridcolor="#F3F4F6",   
                                        showline=True,   
                                        linecolor="#171717",   
                                        linewidth=1,  
                                        ticks="outside",  
                                        tickfont=dict(family="Inter", size=14, color="#171717")  
                                    ),  
                                    plot_bgcolor='white',  
                                    margin=dict(l=10, r=10, t=50, b=10),  
                                    legend=dict(  
                                        font=dict(family="Inter", size=14, color="#171717"),  
                                        orientation="h",  
                                        yanchor="bottom",  
                                        y=1.02,  
                                        xanchor="right",  
                                        x=1  
                                    )  
                                )  
                                # Hide top and right spines  
                                fig_match.update_xaxes(showspikes=True, mirror=False, showline=True, linewidth=1, linecolor="#171717")  
                                fig_match.update_yaxes(showspikes=True, mirror=False, showline=True, linewidth=1, linecolor="#171717")  
                                  
                                st.plotly_chart(fig_match, use_container_width=True)  
                                  
                                # Display matched peaks table  
                                with st.expander("Show Matched Peaks Table"):  
                                    if match_results['match_count'] > 0:  
                                        match_data = {  
                                            "MS2 m/z": [ms2_spec['mz'][i] for i in match_results['ms2_indices']],  
                                            "MS2 Intensity": [ms2_spec['intensity'][i] for i in match_results['ms2_indices']],  
                                            "MGF m/z": [mgf_spec['mz_array'][i] for i in match_results['mgf_indices']],  
                                            "MGF Intensity": [mgf_spec['intensity_array'][i] for i in match_results['mgf_indices']],  
                                            "Difference (Da)": [abs(ms2_spec['mz'][i] - mgf_spec['mz_array'][j])   
                                                              for i, j in zip(match_results['ms2_indices'], match_results['mgf_indices'])]  
                                        }  
                                        match_df = pd.DataFrame(match_data)  
                                        st.dataframe(match_df)  
                                    else:  
                                        st.write("No matching peaks found.")  
                    except Exception as e:  
                        st.error(f"Error parsing MGF file: {str(e)}")  
                        st.info("Please check the format of your MGF file and try again.")  
        except Exception as e:  
            st.error(f"Error loading mzML file: {str(e)}")  
            st.info("Please check the format of your mzML file and try again.")  
else:  
    st.info("Please upload an mzML file to begin analysis.")  
      
    # Show example of what the app can do  
    st.subheader("Features")  
    st.markdown("""  
    - **Total Ion Chromatogram (TIC)**: View the total ion current across retention time  
    - **Extracted Ion Chromatogram (XIC)**: Extract and visualize specific m/z values  
    - **MS2 Fragmentation Analysis**: View MS2 spectra and analyze fragmentation patterns  
    - **MGF Matching**: Match experimental MS2 data with theoretical fragments from MGF files  
    """)  
