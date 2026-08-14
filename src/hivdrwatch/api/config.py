"""API configuration — env-var driven via pydantic-settings.

Anything you might change per environment (rate-limit, DB URL, CORS
origins) lives here. Defaults are safe for local dev.
"""

from __future__ import annotations

from pydantic import Field
from pydantic_settings import BaseSettings, SettingsConfigDict


class APISettings(BaseSettings):
    """Configuration read from env vars prefixed with `HIVDRWATCH_API_`."""

    model_config = SettingsConfigDict(
        env_prefix="HIVDRWATCH_API_",
        env_file=".env",
        env_file_encoding="utf-8",
        extra="ignore",
    )

    # HTTP
    host: str = Field(default="0.0.0.0")
    port: int = Field(default=8080, ge=1, le=65535)
    root_path: str = Field(default="", description="ASGI root_path for reverse proxies")

    # CORS — permissive by default for local dev; lock down in prod.
    cors_origins: list[str] = Field(default_factory=lambda: ["*"])
    cors_allow_credentials: bool = False

    # Rate limit (per IP): N requests per M seconds → e.g. "100/minute"
    rate_limit_default: str = Field(default="100/minute")
    rate_limit_health: str = Field(default="1000/minute")

    # Feature flags
    enable_docs: bool = Field(default=True)


_settings: APISettings | None = None


def get_settings() -> APISettings:
    """Cached settings — reads env once per process."""
    global _settings
    if _settings is None:
        _settings = APISettings()
    return _settings


def reset_settings() -> None:
    """Test-only — force re-read on next call."""
    global _settings
    _settings = None


__all__ = ["APISettings", "get_settings", "reset_settings"]
