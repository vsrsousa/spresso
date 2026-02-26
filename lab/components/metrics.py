"""
Metric cards and statistics components.
"""

import streamlit as st
from ..config import COLORS


def metric_card(label: str, value: str, icon: str = "📊", 
                color: str = None):
    """Display a metric card."""
    if color is None:
        color = COLORS['primary']
    
    st.metric(label=label, value=value)


def stats_row(stats: dict):
    """Display statistics in columns."""
    cols = st.columns(len(stats))
    for col, (label, value) in zip(cols, stats.items()):
        with col:
            st.metric(label=label, value=value)


def info_box(text: str, icon: str = "ℹ️"):
    """Display an information box."""
    st.info(f"{icon} {text}")


def success_box(text: str, icon: str = "✅"):
    """Display a success box."""
    st.success(f"{icon} {text}")


def error_box(text: str, icon: str = "❌"):
    """Display an error box."""
    st.error(f"{icon} {text}")


def warning_box(text: str, icon: str = "⚠️"):
    """Display a warning box."""
    st.warning(f"{icon} {text}")
