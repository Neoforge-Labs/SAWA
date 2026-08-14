import marimo

__generated_with = "0.23.16"
app = marimo.App(width="medium")


@app.cell
def _():
    from Bio import SeqIO
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    import numpy as np
    import polars as pl


    DATA_DIR = "../data/raw/lanl/"
    HIV_FAST_FILE = f"{DATA_DIR}/hiv-db.fasta"
    TSV_FILE = f"{DATA_DIR}/results.tsv"

    lanl_records = []

    with open(HIV_FAST_FILE) as f:
        for header, seq in SimpleFastaParser(f):
            # LANL headers use '.' to separate metadata fields
            parts = header.split('.')
        
            # Guard against unexpected header lengths safely
            if len(parts) >= 5:
                subtype = parts[0]
                country = parts[1]
                year = parts[2]
                name = parts[3]
                accession = parts[4]
            else:
                # Fallback if a header doesn't fit the standard LANL pattern
                subtype, country, year, name, accession = "Unknown", "Unknown", "Unknown", header, "Unknown"
            
            lanl_records.append({
                "subtype": subtype,
                "country": country,
                "year": year,
                "name": name,
                "accession": accession,
                "sequence": seq
            })


    df = pl.DataFrame(lanl_records)
    df

    return df, pl


@app.cell
def _(df):
    country_count = df.group_by(["country", "subtype"]).count().sort("count", descending=True)
    country_count
    return


@app.cell
def _(df):
    subtype_count = df.group_by("subtype").count().sort("count", descending=True)
    subtype_count
    return


@app.cell
def _(df):
    year_count = df.group_by(["year", "subtype"]).count().sort("year", descending=True)
    year_count
    return


@app.cell
def _(df, pl):
    # 1. Clean data fields and handle missing text signals
    cleaned_df = df.with_columns([
        # Standardize all sequence data to uppercase (avoids mixed-case issues)
        pl.col("sequence").str.to_uppercase(),
    
        # Calculate exact sequence length
        pl.col("sequence").str.len_bytes().alias("seq_len"),
    
        # Replace '-' missing year markers with proper Polars nulls, and cast to integer
        pl.when(pl.col("year").str.strip_chars() == "-")
          .then(None)
          .otherwise(pl.col("year"))
          .cast(pl.Int32, strict=False)  # strict=False safely turns unparseable text into null
          .alias("year"),
      
        # Standardize text columns (trim hidden spaces, fix missing flags)
        pl.col("subtype").str.strip_chars().str.to_uppercase(),
        pl.col("country").str.strip_chars().str.to_uppercase(),
    ])

    # 2. Drop duplicates (Prevent counting the exact same record twice)
    # LANL sometimes repeats sequences if they are associated with multiple sub-studies
    cleaned_df = cleaned_df.unique(subset=["accession"])



    print(cleaned_df)
    return (cleaned_df,)


@app.cell
def _(cleaned_df, pl):
    epidemiology_df = cleaned_df.with_columns([
        # Calculate GC Content percentage as a genomic marker
        (((pl.col("sequence").str.count_matches("G") + 
           pl.col("sequence").str.count_matches("C")) / 
          pl.col("seq_len")) * 100).alias("gc_content"),
    
        # Create 5-year cohort windows to track shifting outbreak eras
        ((pl.col("year") // 5) * 5).alias("five_year_window")
    ]).filter(
        pl.col("year").is_not_null() # Remove records missing a timestamp
    )
    epidemiology_df
    return (epidemiology_df,)


@app.cell
def _(epidemiology_df, pl):
    # 1. Calculate historical baseline (average sequences per country/subtype over time)
    baseline_df = (
        epidemiology_df.group_by(["country", "subtype", "year"])
        .len(name="sequence_count")
        .group_by(["country", "subtype"])
        .agg(pl.col("sequence_count").mean().alias("historical_yearly_avg"))
    )

    print(baseline_df)

    # 2. Count sequences in the most recent surveillance data (e.g., your latest year)
    latest_year = epidemiology_df["year"].max()

    latest_counts_df = (
        epidemiology_df.filter(pl.col("year") == latest_year)
        .group_by(["country", "subtype"])
        .len(name="current_year_count")
    )

    print(latest_counts_df)

    # 3. Join and calculate an Outbreak Risk Score
    outbreak_alerts = (
        latest_counts_df.join(baseline_df, on=["country", "subtype"], how="inner")
        .with_columns([
            # Risk Ratio: How many times higher is the current count vs the historical norm?
            (pl.col("current_year_count") / pl.col("historical_yearly_avg")).alias("outbreak_risk_ratio")
        ])
        # Flag high-risk regions where sequences have spiked over 2.5x the baseline
        .filter((pl.col("outbreak_risk_ratio") >= 2.5) & (pl.col("current_year_count") >= 5))
        .sort("outbreak_risk_ratio", descending=True)
    )

    print(f"--- OUTBREAK ALERTS FOR YEAR {latest_year} ---")
    print(outbreak_alerts)

    return


@app.cell
def _(epidemiology_df, pl):
    # Fast K-mer frequency extraction using Polars
    def extract_kmer_features(df, k=3):
        # This slices the sequence string into repeating overlapping chunks of size K
        # E.g., "ATGC" with k=3 becomes ["ATG", "TGC"]
        return df.with_columns([
            pl.col("sequence")
            .str.extract_all(f".{{{k}}}") # Tokenize into chunks of length K
            .alias("kmers")
        ])

    ml_ready_df = extract_kmer_features(epidemiology_df, k=3)

    # View your tokenized dataframe ready for vectorization or frequency tables
    print(ml_ready_df.select(["id", "subtype", "kmers"]).head())

    return


if __name__ == "__main__":
    app.run()
