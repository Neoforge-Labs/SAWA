"""Alembic env — connects Alembic to our SQLModel/SQLAlchemy metadata.

Reads the database URL from `hivdrwatch.storage.session.database_url()`
so `alembic upgrade head` picks up the same URL as the application.
"""

from __future__ import annotations

from logging.config import fileConfig

from alembic import context
from sqlalchemy import engine_from_config, pool

from hivdrwatch.storage.schema import metadata as target_metadata
from hivdrwatch.storage.session import database_url

config = context.config

if config.config_file_name is not None:
    fileConfig(config.config_file_name)

# Override sqlalchemy.url with what session.py resolves to.
config.set_main_option("sqlalchemy.url", database_url())


def run_migrations_offline() -> None:
    """Emit SQL to stdout without a live DB — useful for CI diffs."""
    url = config.get_main_option("sqlalchemy.url")
    context.configure(
        url=url,
        target_metadata=target_metadata,
        literal_binds=True,
        dialect_opts={"paramstyle": "named"},
    )
    with context.begin_transaction():
        context.run_migrations()


def run_migrations_online() -> None:
    """Apply migrations against a live database."""
    connectable = engine_from_config(
        config.get_section(config.config_ini_section, {}),
        prefix="sqlalchemy.",
        poolclass=pool.NullPool,
    )
    with connectable.connect() as connection:
        context.configure(
            connection=connection,
            target_metadata=target_metadata,
        )
        with context.begin_transaction():
            context.run_migrations()


if context.is_offline_mode():
    run_migrations_offline()
else:
    run_migrations_online()
