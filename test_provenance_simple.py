#!/usr/bin/env python3
"""
Simple provenance test - JUST DATABASE INIT
"""

from xespresso.db import DatabaseWorkflow

print("🚀 Testing DatabaseWorkflow initialization...")

# Initialize database workflow
wf = DatabaseWorkflow(db_path='provenance.db')

print("✅ DatabaseWorkflow initialized successfully!")
print(f"   ASE Database: {wf.db_path}")
print(f"   Provenance DB: {wf.provenance.db_path}")

# Test basic provenance operations
cursor = wf.provenance.conn.cursor()
cursor.execute("SELECT COUNT(*) FROM calculations")
count = cursor.fetchone()[0]
print(f"✅ Provenance accessible. Current calculations: {count}")

print("\n🎉 Provenance system is ready!")
print("You can now use DatabaseWorkflow for caching calculations.")