"""initial hivdrwatch schema

Revision ID: 0001_initial
Revises:
Create Date: 2026-08-13 17:00:00.000000
"""

from __future__ import annotations

from collections.abc import Sequence

import sqlalchemy as sa
from alembic import op

revision: str = "0001_initial"
down_revision: str | None = None
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    op.create_table(
        "pipeline_runs",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column(
            "started_at",
            sa.DateTime(timezone=True),
            nullable=False,
            server_default=sa.text("now()"),
        ),
        sa.Column("completed_at", sa.DateTime(timezone=True), nullable=True),
        sa.Column("stage", sa.String(length=32), nullable=False),
        sa.Column("schema_version", sa.String(length=16), nullable=False),
        sa.Column("storage_schema_version", sa.String(length=16), nullable=False),
        sa.Column("git_sha", sa.String(length=40), nullable=True),
        sa.Column("config_hash", sa.String(length=64), nullable=True),
        sa.Column("row_count", sa.Integer(), nullable=True),
        sa.Column("notes", sa.Text(), nullable=True),
    )

    op.create_table(
        "sequences",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("accession", sa.String(length=20), nullable=False),
        sa.Column("schema_version", sa.String(length=16), nullable=False),
        sa.Column("source", sa.String(length=16), nullable=False),
        sa.Column("sequence", sa.Text(), nullable=False),
        sa.Column("length", sa.Integer(), nullable=False),
        sa.Column("gene_region", sa.String(length=8), nullable=False, server_default="pol"),
        sa.Column("country", sa.String(length=2), nullable=True),
        sa.Column("latitude", sa.Float(), nullable=True),
        sa.Column("longitude", sa.Float(), nullable=True),
        sa.Column("collection_year", sa.Integer(), nullable=True),
        sa.Column("subtype", sa.String(length=16), nullable=False, server_default="unknown"),
        sa.Column("subtype_source", sa.String(length=16), nullable=True),
        sa.Column("is_drug_naive", sa.Boolean(), nullable=True),
        sa.Column("risk_group", sa.String(length=32), nullable=True),
        sa.Column("patient_id_hashed", sa.String(length=64), nullable=True),
        sa.Column("pubmed_id", sa.String(length=16), nullable=True),
        sa.Column("study_accession", sa.String(length=32), nullable=True),
        sa.Column("pipeline_run_id", sa.Integer(), nullable=False),
        sa.Column("raw_record_hash", sa.String(length=64), nullable=False),
        sa.Column(
            "ingested_at",
            sa.DateTime(timezone=True),
            nullable=False,
            server_default=sa.text("now()"),
        ),
        sa.UniqueConstraint(
            "schema_version", "accession", name="uq_sequences_schema_accession"
        ),
    )
    op.create_index("ix_sequences_accession", "sequences", ["accession"])
    op.create_index("ix_sequences_country", "sequences", ["country"])
    op.create_index("ix_sequences_collection_year", "sequences", ["collection_year"])
    op.create_index("ix_sequences_subtype", "sequences", ["subtype"])
    op.create_index("ix_sequences_pipeline_run_id", "sequences", ["pipeline_run_id"])

    op.create_table(
        "dr_calls",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("accession", sa.String(length=20), nullable=False),
        sa.Column("drug_class", sa.String(length=8), nullable=False),
        sa.Column("drug", sa.String(length=16), nullable=False),
        sa.Column("level", sa.String(length=16), nullable=False),
        sa.Column("mutations", sa.JSON(), nullable=True),
        sa.Column("algorithm_version", sa.String(length=16), nullable=False),
        sa.Column("pipeline_run_id", sa.Integer(), nullable=False),
        sa.Column(
            "called_at",
            sa.DateTime(timezone=True),
            nullable=False,
            server_default=sa.text("now()"),
        ),
        sa.UniqueConstraint(
            "accession",
            "drug_class",
            "drug",
            "algorithm_version",
            name="uq_dr_calls_seq_drug_algo",
        ),
    )
    op.create_index("ix_dr_calls_accession", "dr_calls", ["accession"])
    op.create_index("ix_dr_calls_drug_class", "dr_calls", ["drug_class"])
    op.create_index("ix_dr_calls_level", "dr_calls", ["level"])
    op.create_index("ix_dr_calls_pipeline_run_id", "dr_calls", ["pipeline_run_id"])

    op.create_table(
        "cluster_memberships",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("accession", sa.String(length=20), nullable=False),
        sa.Column("cluster_id", sa.String(length=32), nullable=False),
        sa.Column("cluster_size", sa.Integer(), nullable=False),
        sa.Column("tn93_threshold", sa.Float(), nullable=False),
        sa.Column("pipeline_run_id", sa.Integer(), nullable=False),
        sa.Column(
            "computed_at",
            sa.DateTime(timezone=True),
            nullable=False,
            server_default=sa.text("now()"),
        ),
        sa.UniqueConstraint(
            "accession",
            "pipeline_run_id",
            name="uq_cluster_memberships_seq_run",
        ),
    )
    op.create_index("ix_cluster_memberships_accession", "cluster_memberships", ["accession"])
    op.create_index("ix_cluster_memberships_cluster_id", "cluster_memberships", ["cluster_id"])
    op.create_index(
        "ix_cluster_memberships_pipeline_run_id", "cluster_memberships", ["pipeline_run_id"]
    )


def downgrade() -> None:
    op.drop_table("cluster_memberships")
    op.drop_table("dr_calls")
    op.drop_table("sequences")
    op.drop_table("pipeline_runs")
