#Author DeepSeek with Anton Smirnov
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import io
from PIL import Image

# Logomaker imports
import logomaker as lm

# Set page config
st.set_page_config(page_title="TCR/BCR Analysis Dashboard", layout="wide")

# ============================================
# HELPER FUNCTIONS
# ============================================

def load_data(uploaded_file):
    """Load CSV data with semicolon delimiter"""
    df = pd.read_csv(uploaded_file, sep=';')
    return df

def calculate_na_statistics(df):
    """Calculate NA statistics for each column"""
    na_stats = pd.DataFrame({
        'Column': df.columns,
        'Missing Values': df.isna().sum(),
        'Missing %': (df.isna().sum() / len(df) * 100).round(2),
        'Present Values': df.count()
    })
    return na_stats

def get_cdr3_statistics(df, filtered_df=None):
    """Calculate comprehensive CDR3 statistics"""
    if filtered_df is None:
        filtered_df = df
    
    stats = {}
    
    # Total number
    stats['total_cdr3_alpha'] = filtered_df['cdr3_alpha'].notna().sum()
    stats['total_cdr3_beta'] = filtered_df['cdr3_beta'].notna().sum()
    stats['total_cdr3'] = stats['total_cdr3_alpha'] + stats['total_cdr3_beta']
    
    # Paired (both chains present)
    stats['paired'] = ((filtered_df['cdr3_alpha'].notna()) & 
                       (filtered_df['cdr3_beta'].notna())).sum()
    
    # By host species
    stats['by_host_species'] = filtered_df.groupby('host_species').agg({
        'cdr3_alpha': lambda x: x.notna().sum(),
        'cdr3_beta': lambda x: x.notna().sum(),
        'id': 'count'
    }).rename(columns={'id': 'total_entries'})
    
    # By chains
    stats['alpha_chains'] = filtered_df[filtered_df['cdr3_alpha'].notna()].shape[0]
    stats['beta_chains'] = filtered_df[filtered_df['cdr3_beta'].notna()].shape[0]
    
    # By epitope
    stats['by_epitope'] = filtered_df.groupby('epitope').agg({
        'cdr3_alpha': lambda x: x.notna().sum(),
        'cdr3_beta': lambda x: x.notna().sum()
    }).fillna(0)
    
    # By database
    stats['by_database'] = filtered_df.groupby('database').agg({
        'cdr3_alpha': lambda x: x.notna().sum(),
        'cdr3_beta': lambda x: x.notna().sum()
    }).fillna(0)
    
    # By MHC allele
    stats['by_mhc'] = filtered_df.groupby(['mhc_alpha', 'mhc_beta']).size().reset_index(name='count')
    stats['by_mhc']['allele'] = stats['by_mhc']['mhc_alpha'].fillna('NA') + '/' + stats['by_mhc']['mhc_beta'].fillna('NA')
    
    return stats

def get_cdr3_per_epitope_data(df, filtered_df=None):
    """Get CDR3 per epitope data without plotting"""
    if filtered_df is None:
        filtered_df = df
    
    # Prepare data for alpha chain
    alpha_data = []
    for idx, row in filtered_df.iterrows():
        if pd.notna(row['cdr3_alpha']) and pd.notna(row['epitope']):
            alpha_data.append({
                'epitope': row['epitope'],
                'cdr3_alpha': row['cdr3_alpha']
            })
    
    # Prepare data for beta chain
    beta_data = []
    for idx, row in filtered_df.iterrows():
        if pd.notna(row['cdr3_beta']) and pd.notna(row['epitope']):
            beta_data.append({
                'epitope': row['epitope'],
                'cdr3_beta': row['cdr3_beta']
            })
    
    alpha_df = pd.DataFrame(alpha_data) if alpha_data else pd.DataFrame()
    beta_df = pd.DataFrame(beta_data) if beta_data else pd.DataFrame()
    
    return alpha_df, beta_df

def get_epitope_statistics(df, filtered_df=None):
    """Calculate epitope statistics"""
    if filtered_df is None:
        filtered_df = df
    
    stats = {}
    
    # Basic counts
    stats['total_epitopes'] = filtered_df['epitope'].nunique()
    stats['total_entries_with_epitope'] = filtered_df['epitope'].notna().sum()
    
    # Top epitopes
    stats['top_epitopes'] = filtered_df['epitope'].value_counts().head(20)
    
    # By epitope source
    stats['by_source'] = filtered_df.groupby('epitope_source').size().reset_index(name='count')
    
    # By species
    stats['by_species'] = filtered_df.groupby('epitope_species').size().reset_index(name='count')
    
    # Epitope length distribution
    epitope_lengths = filtered_df[filtered_df['epitope'].notna()]['epitope'].str.len()
    stats['epitope_lengths'] = epitope_lengths
    
    return stats

def get_vdj_statistics(df, filtered_df=None):
    """Calculate VDJ statistics"""
    if filtered_df is None:
        filtered_df = df
    
    stats = {}
    
    # Count receptors with full info
    stats['full_receptors'] = filtered_df[
        filtered_df['V_alpha'].notna() & 
        filtered_df['V_beta'].notna() & 
        filtered_df['D_beta'].notna() & 
        filtered_df['J_alpha'].notna() & 
        filtered_df['J_beta'].notna() &
        filtered_df['cdr3_alpha'].notna() &
        filtered_df['cdr3_beta'].notna()
    ].shape[0]
    
    # V gene usage
    stats['v_alpha_usage'] = filtered_df['V_alpha'].value_counts().head(20)
    stats['v_beta_usage'] = filtered_df['V_beta'].value_counts().head(20)
    
    # D gene usage (beta only)
    stats['d_beta_usage'] = filtered_df['D_beta'].value_counts().head(20)
    
    # J gene usage
    stats['j_alpha_usage'] = filtered_df['J_alpha'].value_counts().head(20)
    stats['j_beta_usage'] = filtered_df['J_beta'].value_counts().head(20)
    
    # Pairing statistics
    v_pairing = filtered_df.groupby(['V_alpha', 'V_beta']).size().reset_index(name='count')
    stats['v_pairing'] = v_pairing.sort_values('count', ascending=False).head(20)
    
    return stats

def get_mhc_statistics(df, filtered_df=None):
    """Calculate MHC allele statistics"""
    if filtered_df is None:
        filtered_df = df
    
    stats = {}
    
    # By MHC class
    stats['by_class'] = filtered_df.groupby('mhc_class').size().reset_index(name='count')
    
    # Top MHC alleles
    mhc_alleles = filtered_df[filtered_df['mhc_alpha'].notna() & filtered_df['mhc_beta'].notna()]
    if len(mhc_alleles) > 0:
        stats['top_alleles'] = (mhc_alleles['mhc_alpha'] + '/' + mhc_alleles['mhc_beta']).value_counts().head(20)
    else:
        stats['top_alleles'] = pd.Series()
    
    # By MHC alpha
    stats['by_alpha'] = filtered_df['mhc_alpha'].value_counts().head(20)
    
    # By MHC beta
    stats['by_beta'] = filtered_df.query('mhc_beta != "B2M"')['mhc_beta'].value_counts().head(20)
    
    return stats

def create_logomaker_logo(sequences, title="Sequence Logo", max_seqs=500, figsize=(12, 4)):
    """Create professional sequence logo using Logomaker"""
    if len(sequences) == 0:
        return None
    
    # Limit sequences
    sequences = sequences[:max_seqs]
    
    # Get sequence lengths
    lengths = [len(seq) for seq in sequences if isinstance(seq, str)]
    if not lengths:
        return None
    
    max_len = max(lengths)
    
    # Define amino acid order (standard 20 amino acids)
    aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    
    # Initialize count matrix
    counts_matrix = pd.DataFrame(0, index=range(max_len), columns=aa_list)
    
    # Fill matrix
    for seq in sequences:
        if isinstance(seq, str):
            for pos, aa in enumerate(seq[:max_len]):
                if aa in aa_list:
                    counts_matrix.loc[pos, aa] += 1
    
    # Convert to frequencies
    freq_matrix = counts_matrix.div(counts_matrix.sum(axis=1), axis=0)
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    try:
        # Try different color schemes
        logo = lm.Logo(freq_matrix, ax=ax, color_scheme='skylign')
    except:
        try:
            logo = lm.Logo(freq_matrix, ax=ax, color_scheme='chemistry')
        except:
            try:
                logo = lm.Logo(freq_matrix, ax=ax, color_scheme='dmslogo')
            except:
                try:
                    logo = lm.Logo(freq_matrix, ax=ax, color_scheme='classic')
                except:
                    logo = lm.Logo(freq_matrix, ax=ax)
    
    # Customize
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlabel('Position', fontsize=12)
    ax.set_ylabel('Information content (bits)', fontsize=12)
    ax.set_xticks(range(0, max_len, max(1, max_len//10)))
    ax.set_xticklabels(range(1, max_len+1, max(1, max_len//10)))
    
    # Add grid for better readability
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    return fig

# ============================================
# CDR3 PER EPITOPE PLOTS FRAGMENT
# ============================================

@st.fragment
def cdr3_per_epitope_fragment(alpha_df, beta_df, top_n_key="top_n"):
    """Fragment that only updates when the slider changes"""
    
    # Get current top_n from session state or default
    if top_n_key not in st.session_state:
        st.session_state[top_n_key] = 30
    
    # Slider for selecting number of top epitopes
    col1, col2, col3 = st.columns([1, 2, 1])
    with col2:
        top_n = st.slider(
            "Number of top epitopes to display:",
            min_value=5,
            max_value=50,
            value=st.session_state[top_n_key],
            step=5,
            key=top_n_key,
            help="Select how many of the most frequent epitopes to show in the plots"
        )
    
    # Calculate counts with current top_n
    if not alpha_df.empty:
        alpha_counts = alpha_df.groupby('epitope').size().reset_index(name='count')
        alpha_counts = alpha_counts.sort_values('count', ascending=False).head(top_n)
    else:
        alpha_counts = pd.DataFrame()
    
    if not beta_df.empty:
        beta_counts = beta_df.groupby('epitope').size().reset_index(name='count')
        beta_counts = beta_counts.sort_values('count', ascending=False).head(top_n)
    else:
        beta_counts = pd.DataFrame()
    
    # Create two columns for side-by-side plots
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader(f"🔵 Alpha Chain CDR3 per Epitope (Top {top_n})")
        if not alpha_counts.empty:
            fig_alpha = px.bar(
                alpha_counts,
                x='epitope',
                y='count',
                title=f"Number of Alpha CDR3 Sequences per Epitope (Top {top_n})",
                labels={'epitope': 'Epitope', 'count': 'Number of Alpha CDR3 Sequences'},
                color='count',
                color_continuous_scale='Blues',
                text='count'
            )
            fig_alpha.update_traces(textposition='outside', textfont_size=10)
            fig_alpha.update_layout(
                xaxis_tickangle=-45,
                height=500,
                showlegend=False,
                xaxis_title="Epitope",
                yaxis_title="Number of Alpha CDR3 Sequences"
            )
            st.plotly_chart(fig_alpha, use_container_width=True, key="alpha_plot")
            
            # Summary stats for alpha
            col_alpha1, col_alpha2 = st.columns(2)
            with col_alpha1:
                st.metric("Total Alpha CDR3 Sequences", alpha_counts['count'].sum())
            with col_alpha2:
                st.metric("Unique Epitopes (Top {})".format(top_n), len(alpha_counts))
            
            # Download button for alpha data
            csv_alpha = alpha_counts.to_csv(index=False)
            st.download_button(
                label="📥 Download Alpha CDR3 per Epitope Data",
                data=csv_alpha,
                file_name="alpha_cdr3_per_epitope.csv",
                mime="text/csv",
                key="alpha_download"
            )
        else:
            st.info("No alpha chain CDR3 sequences found with associated epitopes")
    
    with col2:
        st.subheader(f"🟢 Beta Chain CDR3 per Epitope (Top {top_n})")
        if not beta_counts.empty:
            fig_beta = px.bar(
                beta_counts,
                x='epitope',
                y='count',
                title=f"Number of Beta CDR3 Sequences per Epitope (Top {top_n})",
                labels={'epitope': 'Epitope', 'count': 'Number of Beta CDR3 Sequences'},
                color='count',
                color_continuous_scale='Greens',
                text='count'
            )
            fig_beta.update_traces(textposition='outside', textfont_size=10)
            fig_beta.update_layout(
                xaxis_tickangle=-45,
                height=500,
                showlegend=False,
                xaxis_title="Epitope",
                yaxis_title="Number of Beta CDR3 Sequences"
            )
            st.plotly_chart(fig_beta, use_container_width=True, key="beta_plot")
            
            # Summary stats for beta
            col_beta1, col_beta2 = st.columns(2)
            with col_beta1:
                st.metric("Total Beta CDR3 Sequences", beta_counts['count'].sum())
            with col_beta2:
                st.metric("Unique Epitopes (Top {})".format(top_n), len(beta_counts))
            
            # Download button for beta data
            csv_beta = beta_counts.to_csv(index=False)
            st.download_button(
                label="📥 Download Beta CDR3 per Epitope Data",
                data=csv_beta,
                file_name="beta_cdr3_per_epitope.csv",
                mime="text/csv",
                key="beta_download"
            )
        else:
            st.info("No beta chain CDR3 sequences found with associated epitopes")
    
    # Optional: Show detailed table for selected epitope
    st.subheader("🔍 View CDR3 Sequences by Epitope and Chain")
    
    # Create options for epitope selection from both chains
    all_epitopes = set()
    if not alpha_counts.empty:
        all_epitopes.update(alpha_counts['epitope'].tolist())
    if not beta_counts.empty:
        all_epitopes.update(beta_counts['epitope'].tolist())
    
    if all_epitopes:
        selected_epitope = st.selectbox(
            "Select an epitope to view its CDR3 sequences:",
            options=sorted(list(all_epitopes)),
            key="epitope_selector"
        )
        
        if selected_epitope:
            col_seq1, col_seq2 = st.columns(2)
            
            # Show alpha sequences
            with col_seq1:
                if not alpha_df.empty:
                    alpha_seqs = alpha_df[alpha_df['epitope'] == selected_epitope]
                    if not alpha_seqs.empty:
                        st.write(f"**Alpha Chain CDR3 sequences for {selected_epitope}:** {len(alpha_seqs)} sequences")
                        st.dataframe(alpha_seqs, key="alpha_seqs_table")
                    else:
                        st.info(f"No alpha CDR3 sequences for {selected_epitope}")
            
            # Show beta sequences
            with col_seq2:
                if not beta_df.empty:
                    beta_seqs = beta_df[beta_df['epitope'] == selected_epitope]
                    if not beta_seqs.empty:
                        st.write(f"**Beta Chain CDR3 sequences for {selected_epitope}:** {len(beta_seqs)} sequences")
                        st.dataframe(beta_seqs, key="beta_seqs_table")
                    else:
                        st.info(f"No beta CDR3 sequences for {selected_epitope}")

# ============================================
# MAIN APP
# ============================================

st.title("🧬 TCR/BCR Sequencing Analysis Dashboard")
st.markdown("Comprehensive analysis of immune receptor repertoire data")

# File upload
uploaded_file = st.sidebar.file_uploader("Upload CSV file", type=['csv'])

if uploaded_file is not None:
    # Load data
    df = load_data(uploaded_file)
    
    st.sidebar.success(f"Loaded {len(df)} entries")
    
    # Show data preview
    with st.expander("Data Preview"):
        st.dataframe(df.head(100))
        st.write(f"Total rows: {len(df)}, Columns: {len(df.columns)}")
    
    # Filters in sidebar
    st.sidebar.markdown("---")
    st.sidebar.header("🔍 Filters")
    
    # Create filtered dataframe based on selections
    filtered_df = df.copy()
    
    # Host species filter
    if 'host_species' in df.columns:
        host_options = ['All'] + list(df['host_species'].dropna().unique())
        selected_host = st.sidebar.selectbox("Host Species", host_options)
        if selected_host != 'All':
            filtered_df = filtered_df[filtered_df['host_species'] == selected_host]
    
    # Database filter
    if 'database' in df.columns:
        db_options = ['All'] + list(df['database'].dropna().unique())
        selected_db = st.sidebar.selectbox("Database", db_options)
        if selected_db != 'All':
            filtered_df = filtered_df[filtered_df['database'] == selected_db]
    
    # Epitope filter
    if 'epitope' in df.columns:
        epitope_options = ['All'] + list(df['epitope'].dropna().unique()[:50])
        selected_epitope = st.sidebar.selectbox("Epitope", epitope_options)
        if selected_epitope != 'All':
            filtered_df = filtered_df[filtered_df['epitope'] == selected_epitope]
    
    # MHC allele filter
    if 'mhc_alpha' in df.columns and 'mhc_beta' in df.columns:
        mhc_options = ['All'] + list((df['mhc_alpha'].fillna('NA') + '/' + df['mhc_beta'].fillna('NA')).unique()[:50])
        selected_mhc = st.sidebar.selectbox("MHC Allele", mhc_options)
        if selected_mhc != 'All':
            mhc_alpha, mhc_beta = selected_mhc.split('/')
            if mhc_alpha != 'NA':
                filtered_df = filtered_df[filtered_df['mhc_alpha'] == mhc_alpha]
            if mhc_beta != 'NA':
                filtered_df = filtered_df[filtered_df['mhc_beta'] == mhc_beta]
    
    # Chain filter for CDR3
    chain_filter = st.sidebar.multiselect(
        "CDR3 Chain",
        options=['alpha', 'beta'],
        default=['alpha', 'beta']
    )
    
    st.sidebar.markdown(f"**Filtered entries:** {len(filtered_df)}")
    
    # ============================================
    # 1. UNIQUE ENTRIES STATISTICS
    # ============================================
    st.header("📊 1. Unique Entries Statistics")
    
    col1, col2 = st.columns(2)
    
    with col1:
        total_unique = len(df.drop_duplicates())
        filtered_unique = len(filtered_df.drop_duplicates())
        st.metric("Total Unique Entries", total_unique)
        st.metric("Filtered Unique Entries", filtered_unique)
    
    with col2:
        if 'host_species' in df.columns:
            unique_by_species = df.groupby('host_species').size().reset_index(name='count')
            fig = px.bar(unique_by_species, x='host_species', y='count', 
                        title="Entries by Host Species")
            st.plotly_chart(fig, use_container_width=True)
    
    # ============================================
    # 2. NA STATISTICS
    # ============================================
    st.header("🔍 2. Missing Values Statistics")
    
    na_stats = calculate_na_statistics(filtered_df)
    fig_na = px.bar(na_stats, x='Column', y='Missing %', 
                    title="Missing Values Percentage by Column")
    st.plotly_chart(fig_na, use_container_width=True)
    st.dataframe(na_stats)
    
    # ============================================
    # 3. CDR3 STATISTICS
    # ============================================
    st.header("🧬 3. CDR3 Statistics")
    
    cdr3_stats = get_cdr3_statistics(df, filtered_df)
    
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Total CDR3 Alpha", cdr3_stats['total_cdr3_alpha'])
    with col2:
        st.metric("Total CDR3 Beta", cdr3_stats['total_cdr3_beta'])
    with col3:
        st.metric("Paired Receptors", cdr3_stats['paired'])
    with col4:
        st.metric("Total CDR3", cdr3_stats['total_cdr3'])
    
    # ============================================
    # CDR3 PER EPITOPE PLOTS (Alpha and Beta) with Fragment
    # ============================================
    st.header("📊 CDR3 Count per Epitope by Chain")
    
    # Get the data once (only recalculates when filters change)
    alpha_df, beta_df = get_cdr3_per_epitope_data(df, filtered_df)
    
    # Call the fragment - only this part will update when slider changes
    cdr3_per_epitope_fragment(alpha_df, beta_df)
    
    # By host species
    st.subheader("CDR3 by Host Species")
    if not cdr3_stats['by_host_species'].empty:
        st.dataframe(cdr3_stats['by_host_species'])
    
    # By epitope
    st.subheader("CDR3 by Epitope (Summary)")
    if not cdr3_stats['by_epitope'].empty:
        st.dataframe(cdr3_stats['by_epitope'])
    
    # By database
    st.subheader("CDR3 by Database")
    if not cdr3_stats['by_database'].empty:
        st.dataframe(cdr3_stats['by_database'])
    
    # By MHC allele
    st.subheader("CDR3 by MHC Allele")
    if not cdr3_stats['by_mhc'].empty:
        st.dataframe(cdr3_stats['by_mhc'].head(20))
    
    # ============================================
    # 4. EPITOPE STATISTICS
    # ============================================
    st.header("🔬 4. Epitope Statistics")
    
    epitope_stats = get_epitope_statistics(df, filtered_df)
    
    col1, col2 = st.columns(2)
    with col1:
        st.metric("Unique Epitopes", epitope_stats['total_epitopes'])
        st.metric("Entries with Epitope", epitope_stats['total_entries_with_epitope'])
    
    with col2:
        if len(epitope_stats['epitope_lengths']) > 0:
            fig_len = px.histogram(epitope_stats['epitope_lengths'], 
                                  title="Epitope Length Distribution",
                                  labels={'value': 'Length', 'count': 'Frequency'})
            st.plotly_chart(fig_len, use_container_width=True)
    
    st.subheader("Top Epitopes")
    st.bar_chart(epitope_stats['top_epitopes'])
    
    # ============================================
    # 5. SEQUENCE LOGOS (with Logomaker)
    # ============================================
    st.header("📈 5. Sequence Logos")
    
    # CDR3 Sequence Logo
    st.subheader("CDR3 Sequence Logo")
    
    # Collect CDR3 sequences based on chain filter
    cdr3_sequences = []
    if 'alpha' in chain_filter:
        alpha_seqs = filtered_df[filtered_df['cdr3_alpha'].notna()]['cdr3_alpha'].tolist()
        cdr3_sequences.extend(alpha_seqs)
    if 'beta' in chain_filter:
        beta_seqs = filtered_df[filtered_df['cdr3_beta'].notna()]['cdr3_beta'].tolist()
        cdr3_sequences.extend(beta_seqs)
    
    if cdr3_sequences:
        fig_logo = create_logomaker_logo(cdr3_sequences, f"CDR3 Sequence Logo ({len(cdr3_sequences)} sequences)")
        if fig_logo:
            st.pyplot(fig_logo)
            plt.close(fig_logo)
        else:
            st.warning("Could not generate CDR3 logo - sequences may be invalid")
    else:
        st.warning("No CDR3 sequences available for logo generation")
    
    # Epitope Sequence Logo
    st.subheader("Epitope Sequence Logo")
    epitope_seqs = filtered_df[filtered_df['epitope'].notna()]['epitope'].tolist()
    
    if epitope_seqs:
        fig_epitope = create_logomaker_logo(epitope_seqs, f"Epitope Sequence Logo ({len(epitope_seqs)} sequences)")
        if fig_epitope:
            st.pyplot(fig_epitope)
            plt.close(fig_epitope)
        else:
            st.warning("Could not generate epitope logo - sequences may be invalid")
    else:
        st.warning("No epitope sequences available for logo generation")
    
    # ============================================
    # 6. VDJ STATISTICS
    # ============================================
    st.header("🧬 6. VDJ Recombination Statistics")
    
    vdj_stats = get_vdj_statistics(df, filtered_df)
    
    st.metric("Receptors with Full VDJ Information", vdj_stats['full_receptors'])
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Top V Alpha Genes")
        if not vdj_stats['v_alpha_usage'].empty:
            fig_v_alpha = px.bar(x=vdj_stats['v_alpha_usage'].values, 
                                 y=vdj_stats['v_alpha_usage'].index,
                                 orientation='h',
                                 title="V Alpha Usage")
            st.plotly_chart(fig_v_alpha, use_container_width=True)
        
        st.subheader("Top D Beta Genes")
        if not vdj_stats['d_beta_usage'].empty:
            fig_d_beta = px.bar(x=vdj_stats['d_beta_usage'].values,
                                y=vdj_stats['d_beta_usage'].index,
                                orientation='h',
                                title="D Beta Usage")
            st.plotly_chart(fig_d_beta, use_container_width=True)
    
    with col2:
        st.subheader("Top V Beta Genes")
        if not vdj_stats['v_beta_usage'].empty:
            fig_v_beta = px.bar(x=vdj_stats['v_beta_usage'].values,
                                y=vdj_stats['v_beta_usage'].index,
                                orientation='h',
                                title="V Beta Usage")
            st.plotly_chart(fig_v_beta, use_container_width=True)
        
        st.subheader("Top J Genes")
        st.write("**J Alpha:**")
        if not vdj_stats['j_alpha_usage'].empty:
            st.bar_chart(vdj_stats['j_alpha_usage'])
        st.write("**J Beta:**")
        if not vdj_stats['j_beta_usage'].empty:
            st.bar_chart(vdj_stats['j_beta_usage'])
    
    st.subheader("Top V Gene Pairings")
    if not vdj_stats['v_pairing'].empty:
        st.dataframe(vdj_stats['v_pairing'])
    
    # ============================================
    # 7. MHC ALLELE STATISTICS
    # ============================================
    st.header("🧬 7. MHC Allele Statistics")
    
    mhc_stats = get_mhc_statistics(df, filtered_df)
    
    col1, col2 = st.columns(2)
    
    with col1:
        if not mhc_stats['by_class'].empty:
            fig_class = px.pie(mhc_stats['by_class'], values='count', names='mhc_class',
                              title="MHC Class Distribution")
            st.plotly_chart(fig_class, use_container_width=True)
        
        st.subheader("Top MHC Alpha Alleles")
        if not mhc_stats['by_alpha'].empty:
            st.bar_chart(mhc_stats['by_alpha'].head(15))
    
    with col2:
        st.subheader("Top MHC Beta Alleles")
        if not mhc_stats['by_beta'].empty:
            st.bar_chart(mhc_stats['by_beta'].head(15))
    
    st.subheader("Top Combined MHC Alleles")
    if not mhc_stats['top_alleles'].empty:
        st.bar_chart(mhc_stats['top_alleles'].head(20))
    
    # Download filtered data
    st.markdown("---")
    csv = filtered_df.to_csv(index=False, sep=';')
    st.download_button(
        label="📥 Download Filtered Data",
        data=csv,
        file_name="filtered_data.csv",
        mime="text/csv"
    )

else:
    st.info("👈 Please upload a CSV file to begin analysis")
    st.markdown("""
    ### Expected CSV format with semicolon delimiter:
    - id
    - D_alpha, D_beta
    - J_alpha, J_beta
    - V_alpha, V_beta
    - cdr3_alpha, cdr3_beta
    - epitope
    - mhc_alpha, mhc_beta
    - mhc_class
    - host_species
    - epitope_species
    - epitope_source
    - database
    """)