#!/usr/bin/env python3

"""
MiGEx Report Generator
Generates comprehensive PDF reports from MiGEx pipeline results
"""

import sys
import os
import json
import pandas as pd
from datetime import datetime
from pathlib import Path
from io import BytesIO
import zipfile

from reportlab.lib.pagesizes import letter, A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, PageBreak, Image
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT, TA_JUSTIFY


def extract_fastqc_data(zip_path):
    """Extract quality metrics from FastQC data file within a zip archive"""
    qc_metrics = {
        'sequences': 'N/A',
        'poor_quality': 'N/A',
        'sequence_length': 'N/A',
        'gc_percent': 'N/A',
        'summary': {}
    }
    
    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            # Find fastqc_data.txt in the zip
            for file_info in zip_ref.filelist:
                if file_info.filename.endswith('fastqc_data.txt'):
                    with zip_ref.open(file_info.filename) as f:
                        content = f.read().decode('utf-8')
                        for line in content.split('\n'):
                            if line.startswith('Total Sequences'):
                                qc_metrics['sequences'] = line.split('\t')[1].strip()
                            elif line.startswith('Sequence length'):
                                qc_metrics['sequence_length'] = line.split('\t')[1].strip()
                            elif line.startswith('%GC'):
                                qc_metrics['gc_percent'] = line.split('\t')[1].strip()
                            elif line.startswith('Sequences flagged as poor quality'):
                                qc_metrics['poor_quality'] = line.split('\t')[1].strip()
                            elif line.startswith('>>'):
                                # Parse summary section
                                parts = line[2:].split('\t')
                                if len(parts) >= 2:
                                    category = parts[0].strip()
                                    status = parts[1].strip().lower()
                                    qc_metrics['summary'][category] = status
                    break
    except Exception as e:
        pass
    
    return qc_metrics


def extract_qc_metrics(results_dir, sample):
    """Extract QC metrics from FastQC result files"""
    qc_data = {}
    
    raw_qc_dir = os.path.join(results_dir, "1-QC/FASTQC_raw")
    filtered_qc_dir = os.path.join(results_dir, "1-QC/FASTQC_filtered")
    
    # Extract metrics for forward and reverse reads
    for read_num in ['1', '2']:
        read_label = 'Forward' if read_num == '1' else 'Reverse'
        
        # Find raw fastqc zip file
        raw_zip = None
        if os.path.exists(raw_qc_dir):
            for f in os.listdir(raw_qc_dir):
                if f.endswith('.zip') and f"_{read_num}_" in f:
                    raw_zip = os.path.join(raw_qc_dir, f)
                    break
        
        # Find filtered fastqc zip file
        filtered_zip = None
        if os.path.exists(filtered_qc_dir):
            for f in os.listdir(filtered_qc_dir):
                if f.endswith('.zip') and f"_{read_num}.filtered" in f:
                    filtered_zip = os.path.join(filtered_qc_dir, f)
                    break
        
        # Extract metrics
        raw_metrics = extract_fastqc_data(raw_zip) if raw_zip else {}
        filtered_metrics = extract_fastqc_data(filtered_zip) if filtered_zip else {}
        
        # Store data for comparison
        qc_data[f'{read_label} Read - Raw'] = {
            'Number of Sequences': raw_metrics.get('sequences', 'N/A'),
            'Sequence Length': raw_metrics.get('sequence_length', 'N/A'),
            '%GC': raw_metrics.get('gc_percent', 'N/A'),
            'Poor Quality': raw_metrics.get('poor_quality', 'N/A'),
            'summary': raw_metrics.get('summary', {})
        }
        
        qc_data[f'{read_label} Read - Filtered'] = {
            'Number of Sequences': filtered_metrics.get('sequences', 'N/A'),
            'Sequence Length': filtered_metrics.get('sequence_length', 'N/A'),
            '%GC': filtered_metrics.get('gc_percent', 'N/A'),
            'Poor Quality': filtered_metrics.get('poor_quality', 'N/A'),
            'summary': filtered_metrics.get('summary', {})
        }
    
    return qc_data


def extract_assembly_metrics(results_dir, sample):
    """Extract assembly metrics from QUAST report"""
    assembly_data = {}
    quast_report = os.path.join(results_dir, f"2-Assembly/quast/{sample}/report.tsv")
    
    if os.path.exists(quast_report):
        df = pd.read_csv(quast_report, sep='\t', header=None)
        for idx, row in df.iterrows():
            metric = row[0]
            value = row[1]
            assembly_data[metric] = str(value)
    
    return assembly_data


def extract_coverage_metrics(results_dir, sample):
    """Extract coverage metrics from QualiMap"""
    coverage_data = {}
    qualimap_dir = os.path.join(results_dir, f"3-Analysis/coverage/{sample}_qualimap")
    
    if os.path.exists(qualimap_dir):
        # Try to read genome_results.txt from QualiMap
        results_file = os.path.join(qualimap_dir, "genome_results.txt")
        if os.path.exists(results_file):
            with open(results_file, 'r') as f:
                for line in f:
                    if '=' in line:
                        key, value = line.strip().split('=', 1)
                        coverage_data[key.strip()] = value.strip()
    
    return coverage_data


def extract_amr_genes(results_dir, sample):
    """Extract AMR genes from AMRfinder"""
    amr_data = []
    amr_file = os.path.join(results_dir, f"3-Analysis/amrfinder/{sample}_amrfinder.tsv")
    
    if os.path.exists(amr_file):
        df = pd.read_csv(amr_file, sep='\t')
        if len(df) > 0:
            # Select relevant columns
            cols_to_show = ['Gene symbol', 'Protein name', 'Class', 'Subclass', 'Coverage', 'Identity']
            available_cols = [c for c in cols_to_show if c in df.columns]
            amr_data = df[available_cols].head(20).values.tolist()  # Limit to 20 rows
    
    return amr_data


def extract_virulence_factors(results_dir, sample):
    """Extract virulence factors from ABRicate"""
    vf_data = []
    vf_file = os.path.join(results_dir, f"3-Analysis/abricate/{sample}_virulence.tsv")
    
    if os.path.exists(vf_file):
        df = pd.read_csv(vf_file, sep='\t')
        if len(df) > 0:
            # Extract relevant columns
            cols_to_show = ['GENE', 'IDENTITY', 'COVERAGE', 'DATABASE']
            available_cols = [c for c in cols_to_show if c in df.columns]
            if available_cols:
                vf_data = df[available_cols].head(15).values.tolist()
    
    return vf_data


def extract_plasmids(results_dir, sample):
    """Extract plasmid info from ABRicate"""
    plasmid_data = []
    plasmid_file = os.path.join(results_dir, f"3-Analysis/abricate/{sample}_plasmids.tsv")
    
    if os.path.exists(plasmid_file):
        df = pd.read_csv(plasmid_file, sep='\t')
        if len(df) > 0:
            cols_to_show = ['GENE', 'IDENTITY', 'COVERAGE', 'DATABASE']
            available_cols = [c for c in cols_to_show if c in df.columns]
            if available_cols:
                plasmid_data = df[available_cols].values.tolist()
    
    return plasmid_data


def extract_mlst_type(results_dir, sample):
    """Extract MLST type"""
    mlst_data = {}
    mlst_file = os.path.join(results_dir, f"3-Analysis/mlst/{sample}_mlst.tsv")
    
    if os.path.exists(mlst_file):
        df = pd.read_csv(mlst_file, sep='\t')
        if len(df) > 0:
            row = df.iloc[0]
            mlst_data['Scheme'] = str(row.get('scheme', 'N/A'))
            mlst_data['ST (Sequence Type)'] = str(row.get('ST', 'N/A'))
            # Get allele info
            allele_cols = [c for c in df.columns if c not in ['FILE', 'scheme', 'ST', 'aLocus']]
            if allele_cols:
                mlst_data['Alleles'] = ' | '.join([f"{c}:{row.get(c, 'N/A')}" for c in allele_cols[:5]])
    
    return mlst_data


def extract_taxonomy(results_dir, sample):
    """Extract taxonomy from Mash"""
    tax_data = []
    mash_file = os.path.join(results_dir, f"3-Analysis/mash/{sample}_mash_distances.txt")
    
    if os.path.exists(mash_file):
        with open(mash_file, 'r') as f:
            lines = f.readlines()
            for line in lines[:10]:  # Top 10 matches
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    tax_data.append([parts[0].split('/')[-1].replace('.gz', ''), 
                                    f"{float(parts[2]):.6f}", 
                                    f"{parts[3] if len(parts) > 3 else 'N/A'}"])
    
    return tax_data


def extract_bakta_stats(results_dir, sample):
    """Extract gene annotation stats from Bakta"""
    stats = {}
    bakta_json = os.path.join(results_dir, f"3-Analysis/bakta/{sample}/{sample}.json")
    
    if os.path.exists(bakta_json):
        with open(bakta_json, 'r') as f:
            data = json.load(f)
            stats['Total genes'] = str(data.get('statistics', {}).get('gene', 'N/A'))
            stats['CDSs'] = str(data.get('statistics', {}).get('cds', 'N/A'))
            stats['rRNA'] = str(data.get('statistics', {}).get('rrna', 'N/A'))
            stats['tRNA'] = str(data.get('statistics', {}).get('trna', 'N/A'))
    
    return stats


def extract_qualimap_figure(results_dir, sample):
    """Extract QualiMap HTML report path"""
    qualimap_html = os.path.join(results_dir, f"3-Analysis/coverage/{sample}_qualimap/qualimapReport.html")
    return qualimap_html if os.path.exists(qualimap_html) else None


def create_pdf_report(results_dir, sample, output_dir):
    """Create comprehensive PDF report"""
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # PDF output path
    pdf_path = os.path.join(output_dir, f"{sample}_report.pdf")
    
    # Create PDF
    doc = SimpleDocTemplate(pdf_path, pagesize=letter,
                           topMargin=0.5*inch, bottomMargin=0.5*inch,
                           leftMargin=0.75*inch, rightMargin=0.75*inch)
    
    # Container for PDF elements
    elements = []
    styles = getSampleStyleSheet()
    
    # Custom styles
    title_style = ParagraphStyle(
        'CustomTitle',
        parent=styles['Heading1'],
        fontSize=24,
        textColor=colors.HexColor('#1f4788'),
        spaceAfter=12,
        alignment=TA_CENTER,
        fontName='Helvetica-Bold'
    )
    
    subtitle_style = ParagraphStyle(
        'CustomSubtitle',
        parent=styles['Heading3'],
        fontSize=14,
        textColor=colors.HexColor('#2e5c99'),
        spaceAfter=12,
        alignment=TA_CENTER,
        fontName='Helvetica'
    )
    
    heading_style = ParagraphStyle(
        'CustomHeading',
        parent=styles['Heading2'],
        fontSize=14,
        textColor=colors.HexColor('#2e5c99'),
        spaceAfter=8,
        spaceBefore=12,
        fontName='Helvetica-Bold'
    )
    
    # ===== COVER PAGE =====
    elements.append(Spacer(1, 1.5*inch))
    elements.append(Paragraph("MiGEx Report", title_style))
    elements.append(Paragraph("Microbial Genome Explorer", subtitle_style))
    elements.append(Spacer(1, 0.5*inch))
    
    cover_data = [
        ['Sample ID:', sample],
        ['Report Date:', datetime.now().strftime("%Y-%m-%d %H:%M:%S")],
        ['Pipeline Version:', '0.1.0']
    ]
    cover_table = Table(cover_data, colWidths=[2*inch, 3*inch])
    cover_table.setStyle(TableStyle([
        ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
        ('FONTNAME', (0, 0), (0, -1), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, -1), 10),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
    ]))
    elements.append(cover_table)
    elements.append(PageBreak())
    
    # ===== QC SUMMARY =====
    elements.append(Paragraph("1. Quality Control Summary", heading_style))
    qc_metrics = extract_qc_metrics(results_dir, sample)
    
    if qc_metrics:
        # Compact heading style
        compact_heading = ParagraphStyle(
            'CompactHeading',
            parent=styles['Heading4'],
            fontSize=12,
            spaceAfter=6,
            spaceBefore=6
        )
        
        # Forward reads comparison
        elements.append(Paragraph("Reads Comparison (R1/R2)", compact_heading))
        forward_data = [['Read', 'Condition', '# Seq', 'Length', '%GC', 'Poor Quality']]
        if 'Forward Read - Raw' in qc_metrics:
            raw_row = qc_metrics['Forward Read - Raw']
            forward_data.append(['R1', 'Raw', raw_row.get('Number of Sequences', 'N/A'), 
                               raw_row.get('Sequence Length', 'N/A'), 
                               raw_row.get('%GC', 'N/A'), 
                               raw_row.get('Poor Quality', 'N/A')])
        if 'Forward Read - Filtered' in qc_metrics:
            filtered_row = qc_metrics['Forward Read - Filtered']
            forward_data.append(['R1', 'Filtered', filtered_row.get('Number of Sequences', 'N/A'), 
                               filtered_row.get('Sequence Length', 'N/A'), 
                               filtered_row.get('%GC', 'N/A'), 
                               filtered_row.get('Poor Quality', 'N/A')])
        if 'Reverse Read - Raw' in qc_metrics:
            raw_row = qc_metrics['Reverse Read - Raw']
            forward_data.append(['R2', 'Raw', raw_row.get('Number of Sequences', 'N/A'), 
                               raw_row.get('Sequence Length', 'N/A'), 
                               raw_row.get('%GC', 'N/A'), 
                               raw_row.get('Poor Quality', 'N/A')])
        if 'Reverse Read - Filtered' in qc_metrics:
            filtered_row = qc_metrics['Reverse Read - Filtered']
            forward_data.append(['R2', 'Filtered', filtered_row.get('Number of Sequences', 'N/A'), 
                               filtered_row.get('Sequence Length', 'N/A'), 
                               filtered_row.get('%GC', 'N/A'), 
                               filtered_row.get('Poor Quality', 'N/A')])
        
        forward_table = Table(forward_data, colWidths=[0.5*inch, 0.8*inch, 1.2*inch, 1*inch, 0.7*inch])
        forward_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#2e5c99')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 7),
            ('TOPPADDING', (0, 0), (-1, -1), 3),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 3),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f8f8f8')])
        ]))
        elements.append(forward_table)
        elements.append(Spacer(1, 0.08*inch))
        
        # FastQC Summary Flags
        elements.append(Paragraph("Quality Flags", compact_heading))
        
        # Prepare summary data
        def get_status_color(status):
            """Return color based on status"""
            status_lower = status.lower()
            if status_lower == 'pass':
                return colors.HexColor('#90EE90')  # Light green
            elif status_lower == 'warn':
                return colors.HexColor('#FFD700')  # Gold/Orange
            elif status_lower == 'fail':
                return colors.HexColor('#FF6B6B')  # Light red
            return colors.white
        
        # Combine all summaries into one table
        all_summaries = {
            'R1 Raw': qc_metrics.get('Forward Read - Raw', {}).get('summary', {}),
            'R1 Filtered': qc_metrics.get('Forward Read - Filtered', {}).get('summary', {}),
            'R2 Raw': qc_metrics.get('Reverse Read - Raw', {}).get('summary', {}),
            'R2 Filtered': qc_metrics.get('Reverse Read - Filtered', {}).get('summary', {})
        }
        
        # Get all categories
        all_categories = set()
        for summary in all_summaries.values():
            all_categories.update(summary.keys())
        
        summary_data = [['Category', 'R1 Raw', 'R1 Filt', 'R2 Raw', 'R2 Filt']]
        for category in sorted(all_categories):
            row = [category]
            for key in ['R1 Raw', 'R1 Filtered', 'R2 Raw', 'R2 Filtered']:
                row.append(all_summaries[key].get(category, 'N/A'))
            summary_data.append(row)
        
        if len(summary_data) > 1:
            summary_table = Table(summary_data, colWidths=[2.2*inch, 0.7*inch, 0.7*inch, 0.7*inch, 0.7*inch])
            table_style = [
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#2e5c99')),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, -1), 6),
                ('TOPPADDING', (0, 0), (-1, -1), 2),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 2),
                ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
            ]
            
            # Add colored backgrounds for each status
            for row_idx in range(1, len(summary_data)):
                for col_idx in range(1, 5):
                    status = summary_data[row_idx][col_idx]
                    table_style.append(('BACKGROUND', (col_idx, row_idx), (col_idx, row_idx), get_status_color(status)))
            
            summary_table.setStyle(TableStyle(table_style))
            elements.append(summary_table)
    
    elements.append(PageBreak())
    
    # ===== ASSEMBLY METRICS =====
    elements.append(Paragraph("2. Assembly Statistics", heading_style))
    assembly_metrics = extract_assembly_metrics(results_dir, sample)
    
    if assembly_metrics:
        assembly_data = [['Metric', 'Value']]
        for k, v in list(assembly_metrics.items())[:15]:  # Top 15 metrics
            assembly_data.append([k, v])
        
        assembly_table = Table(assembly_data, colWidths=[3*inch, 2*inch])
        assembly_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#2e5c99')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('ALIGN', (1, 1), (1, -1), 'RIGHT'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 9),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
            ('GRID', (0, 0), (-1, -1), 1, colors.grey),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f0f0f0')])
        ]))
        elements.append(assembly_table)
    
    elements.append(PageBreak())
    
    # ===== COVERAGE ANALYSIS =====
    elements.append(Paragraph("3. Coverage Analysis", heading_style))
    coverage_metrics = extract_coverage_metrics(results_dir, sample)
    
    if coverage_metrics:
        coverage_data = [['Metric', 'Value']]
        for k, v in list(coverage_metrics.items())[:10]:
            coverage_data.append([k[:50], str(v)[:50]])  # Truncate long values
        
        coverage_table = Table(coverage_data, colWidths=[2.5*inch, 2.5*inch])
        coverage_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#2e5c99')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 8),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
            ('GRID', (0, 0), (-1, -1), 1, colors.grey),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f0f0f0')])
        ]))
        elements.append(coverage_table)
    
    elements.append(PageBreak())
    
    # ===== GENOMIC FEATURES =====
    elements.append(Paragraph("4. Genomic Features", heading_style))
    bakta_stats = extract_bakta_stats(results_dir, sample)
    
    if bakta_stats:
        stats_data = [['Feature', 'Count']]
        for k, v in bakta_stats.items():
            stats_data.append([k, v])
        
        stats_table = Table(stats_data, colWidths=[2.5*inch, 2*inch])
        stats_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#2e5c99')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('ALIGN', (1, 1), (1, -1), 'RIGHT'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 9),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
            ('GRID', (0, 0), (-1, -1), 1, colors.grey),
        ]))
        elements.append(stats_table)
    
    elements.append(PageBreak())
    
    # ===== ANTIMICROBIAL RESISTANCE =====
    elements.append(Paragraph("5. Antimicrobial Resistance Genes", heading_style))
    amr_genes = extract_amr_genes(results_dir, sample)
    
    if amr_genes:
        amr_data = [['Gene', 'Protein', 'Class', 'Coverage', 'Identity']]
        for gene in amr_genes[:10]:  # Top 10
            amr_data.append([str(g)[:20] for g in gene[:5]])
        
        amr_table = Table(amr_data, colWidths=[1.2*inch, 1.5*inch, 1.2*inch, 0.8*inch, 0.8*inch])
        amr_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#c1504e')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 8),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
            ('GRID', (0, 0), (-1, -1), 1, colors.grey),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#ffe6e6')])
        ]))
        elements.append(amr_table)
    else:
        elements.append(Paragraph("No antimicrobial resistance genes detected.", styles['Normal']))
    
    elements.append(PageBreak())
    
    # ===== VIRULENCE FACTORS =====
    elements.append(Paragraph("6. Virulence Factors", heading_style))
    vf_factors = extract_virulence_factors(results_dir, sample)
    
    if vf_factors:
        vf_data = [['Gene', 'Identity', 'Coverage', 'Database']]
        for vf in vf_factors[:10]:
            vf_data.append([str(v)[:25] for v in vf])
        
        vf_table = Table(vf_data, colWidths=[2*inch, 1*inch, 1*inch, 1.5*inch])
        vf_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#f7941d')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 8),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
            ('GRID', (0, 0), (-1, -1), 1, colors.grey),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#fff5eb')])
        ]))
        elements.append(vf_table)
    else:
        elements.append(Paragraph("No virulence factors detected.", styles['Normal']))
    
    elements.append(PageBreak())
    
    # ===== PLASMIDS =====
    elements.append(Paragraph("7. Plasmid Elements", heading_style))
    plasmids = extract_plasmids(results_dir, sample)
    
    if plasmids:
        plasmid_data = [['Element', 'Identity', 'Coverage', 'Database']]
        for p in plasmids[:10]:
            plasmid_data.append([str(x)[:25] for x in p])
        
        plasmid_table = Table(plasmid_data, colWidths=[2*inch, 1*inch, 1*inch, 1.5*inch])
        plasmid_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#7030a0')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 8),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
            ('GRID', (0, 0), (-1, -1), 1, colors.grey),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f3e6ff')])
        ]))
        elements.append(plasmid_table)
    else:
        elements.append(Paragraph("No plasmid elements detected.", styles['Normal']))
    
    elements.append(PageBreak())
    
    # ===== STRAIN IDENTIFICATION =====
    elements.append(Paragraph("8. Strain Identification", heading_style))
    
    # MLST
    elements.append(Paragraph("MLST (Multi-Locus Sequence Type)", styles['Heading4']))
    mlst_type = extract_mlst_type(results_dir, sample)
    
    if mlst_type:
        mlst_data = [['Metric', 'Value']]
        for k, v in mlst_type.items():
            mlst_data.append([k, str(v)[:50]])
        
        mlst_table = Table(mlst_data, colWidths=[2*inch, 3*inch])
        mlst_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#2e5c99')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 9),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
            ('GRID', (0, 0), (-1, -1), 1, colors.grey),
        ]))
        elements.append(mlst_table)
    else:
        elements.append(Paragraph("MLST type could not be determined.", styles['Normal']))
    
    elements.append(Spacer(1, 0.3*inch))
    
    # Taxonomic confirmation
    elements.append(Paragraph("Taxonomic Confirmation (Mash)", styles['Heading4']))
    taxonomy = extract_taxonomy(results_dir, sample)
    
    if taxonomy:
        tax_data = [['Reference Genome', 'Distance', 'P-value']]
        for tax in taxonomy[:5]:
            tax_data.append(tax)
        
        tax_table = Table(tax_data, colWidths=[3*inch, 1.2*inch, 1.3*inch])
        tax_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#2e5c99')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 8),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
            ('GRID', (0, 0), (-1, -1), 1, colors.grey),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f0f0f0')])
        ]))
        elements.append(tax_table)
    else:
        elements.append(Paragraph("Taxonomic information not available.", styles['Normal']))
    
    elements.append(PageBreak())
    
    # ===== FOOTER =====
    elements.append(Spacer(1, 1*inch))
    footer_text = f"Report generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}<br/>MiGEx Pipeline v0.1.0"
    elements.append(Paragraph(footer_text, ParagraphStyle(
        'Footer',
        parent=styles['Normal'],
        fontSize=8,
        textColor=colors.grey,
        alignment=TA_CENTER
    )))
    
    # Build PDF
    try:
        doc.build(elements)
        print(f"✓ Report successfully generated: {pdf_path}")
        return True
    except Exception as e:
        print(f"✗ Error generating report: {e}")
        return False


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 generate_report.py <results_directory> [sample_name] [output_directory]")
        print("Example: python3 generate_report.py /path/to/results SRR5006289 /path/to/reports")
        sys.exit(1)
    
    results_dir = sys.argv[1]
    sample = sys.argv[2] if len(sys.argv) > 2 else None
    output_dir = sys.argv[3] if len(sys.argv) > 3 else os.path.join(results_dir, "4-Reports")
    
    if not sample:
        print("Error: sample name required")
        sys.exit(1)
    
    if not os.path.exists(results_dir):
        print(f"Error: results directory not found: {results_dir}")
        sys.exit(1)
    
    success = create_pdf_report(results_dir, sample, output_dir)
    sys.exit(0 if success else 1)
