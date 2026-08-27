"""Interfaccia web locale. Il backend e' un CLIENT delle stesse funzioni della CLI."""
from zefiro.web.app import create_app  # noqa: F401
from zefiro.web.jobs import JobRegistry  # noqa: F401
