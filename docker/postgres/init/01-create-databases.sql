-- Postgres init: create separate databases for Prefect and MLflow.
--
-- Runs once at container first boot from
-- /docker-entrypoint-initdb.d/. If POSTGRES_USER differs from the
-- default `hivdrwatch`, adjust the OWNER clauses below.

CREATE DATABASE prefect;
CREATE DATABASE mlflow;

-- No separate app users for now — the single POSTGRES_USER owns both
-- databases. Split into per-service users if you ever expose Postgres
-- beyond localhost.
