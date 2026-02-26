"""
Caching utilities for XESPRESSO Lab.
"""

import streamlit as st
from functools import wraps
from .config import CACHE_TTL


def clear_all_caches():
    """Clear all Streamlit caches."""
    st.cache_data.clear()
    st.cache_resource.clear()


def cached_query(ttl: int = CACHE_TTL):
    """Decorator for caching query results."""
    def decorator(func):
        @wraps(func)
        @st.cache_data(ttl=ttl)
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)
        return wrapper
    return decorator
