import marimo

__generated_with = "0.23.16"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import polars as pl
    import numpy as np
    from dowhy import CausalModel
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.metrics import classification_report

    # Marimo Interactive UI widgets
    region_selector = mo.ui.dropdown(
        options=["Kilifi", "Turkana", "Kwale", "All"], 
        value="All", 
        label="Select Clinic Region:"
    )
    region_selector

    return (
        CausalModel,
        RandomForestClassifier,
        classification_report,
        mo,
        np,
        pl,
        region_selector,
    )


@app.cell
def _(mo, np, pl, region_selector):
    # This cell runs automatically whenever 'region_selector.value' changes
    selected_region = region_selector.value

    # 1. Ingesting raw health data using Polars
    # (Simulated data reflecting low-resource tropical disease dynamics)
    np.random.seed(42)
    n_patients = 2000

    df_raw = pl.DataFrame({
        "Malnutrition": np.random.binomial(1, 0.45 if selected_region == "Turkana" else 0.30, n_patients),
        "Chronic_Helminths": np.random.binomial(1, 0.35 if selected_region == "Kilifi" else 0.20, n_patients),
        "Target_Tropical_Disease": np.random.binomial(1, 0.25, n_patients),
        "Hb_Count": np.random.uniform(5.0, 16.0, n_patients),
        "WBC_Count": np.random.uniform(2.0, 20.0, n_patients),
        "Patient_Survival": np.random.binomial(1, 0.85, n_patients)
    })

    mo.md(f"Loaded **{n_patients}** patient records for region: **{selected_region}**.")

    return (df_raw,)


@app.cell
def _(CausalModel, df_raw, mo, pl):
    # Define assumptions using the DOT graph language we explored
    fbc_dag_dot = """
    digraph {
        Malnutrition -> Low_Hemoglobin; Malnutrition -> Patient_Survival;
        Chronic_Helminths -> Low_Hemoglobin; Chronic_Helminths -> WBC_Abnormality; Chronic_Helminths -> Patient_Survival;
        Target_Tropical_Disease -> Low_Hemoglobin; Target_Tropical_Disease -> WBC_Abnormality; Target_Tropical_Disease -> Patient_Survival;
        Low_Hemoglobin -> Patient_Survival; WBC_Abnormality -> Patient_Survival;
    }
    """

    # Map to Pandas inside the cell block for DoWhy's engine compatibility
    df_pd_bridge = df_raw.with_columns([
        (pl.col("Hb_Count") < 11.0).cast(pl.Float64).alias("Low_Hemoglobin"),
        (pl.col("WBC_Count") > 11.0).cast(pl.Float64).alias("WBC_Abnormality")
    ]).to_pandas()

    causal_model = CausalModel(
        data=df_pd_bridge,
        treatment="Low_Hemoglobin",
        outcome="Patient_Survival",
        graph=fbc_dag_dot.replace("\n", " ")
    )

    identified_estimand = causal_model.identify_effect(proceed_when_unidentifiable=True)

    # Format the output beautifully inside Marimo using Markdown
    mo.vstack([
        mo.md("### 🛡️ Causal Guardrail Identifications"),
        mo.plain(identified_estimand)
    ])

    return


@app.cell
def _(df_raw, mo, pl):
    # Use Polars expressions to binarize continuous laboratory blood metrics
    df_features = df_raw.with_columns([
        pl.when(pl.col("Hb_Count") < 11.0).then(1.0).otherwise(0.0).alias("Low_Hemoglobin"),
        pl.when(pl.col("WBC_Count") > 11.0).then(1.0).otherwise(0.0).alias("WBC_Abnormality")
    ]).select([
        # The exact backdoor features identified by DoWhy + our blood panel indicators
        "Low_Hemoglobin", 
        "WBC_Abnormality", 
        "Malnutrition", 
        "Chronic_Helminths", 
        "Target_Tropical_Disease", # Target class label for disease classification
    ])

    # Display a clean, interactive preview table directly in the Marimo UI
    mo.ui.table(df_features.head(5))

    return (df_features,)


@app.cell
def _(RandomForestClassifier, classification_report, df_features, mo):
    # Prepare matrix inputs for Scikit-Learn
    X = df_features.select(["Low_Hemoglobin", "WBC_Abnormality", "Malnutrition", "Chronic_Helminths"]).to_numpy()
    y = df_features.select("Target_Tropical_Disease").to_numpy().ravel()

    # FIX: Added class_weight="balanced" so the tree penalizes missing the minority sick class
    clf = RandomForestClassifier(max_depth=5, random_state=42, class_weight="balanced")
    clf.fit(X, y)

    # Evaluate performance
    predictions = clf.predict(X)
    report = classification_report(y, predictions, output_dict=True)

    # Generate a clean, reactive dashboard view using your precise dictionary keys
    mo.hstack([
        mo.stat(
            value=f"{report['accuracy']:.2%}", 
            label="Overall Diagnostic Accuracy"
        ),
        mo.stat(
            value=f"{report['1']['precision']:.2%}", 
            label="Tropical Disease Precision Rate"
        ),
        mo.stat(
            value=f"{report['1']['recall']:.2%}", 
            label="Sensitivity / Recall (Catch Rate)"
        )
    ], justify="space-around")

    return


if __name__ == "__main__":
    app.run()
