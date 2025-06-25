from setuptools import setup, find_packages

setup(
    name="quickortep",
    version="0.1.0",
    description="A lightweight, fast molecular visualization tool with interactive bond editing",
    packages=find_packages(),
    py_modules=[
        "quickORTEP", "bond_edit_tracker", "bond_manager", "bonds", "config", 
        "elements_table", "event_handler", "export", "geometry_utils", 
        "graph_viewer", "help_text", "hud", "message_panel", "message_service", 
        "molecule", "molecule_nci", "normal_modes", "ortep_molecule", 
        "ortep_renderer", "ortep_viewer", "selection_manager", "test_shapes", 
        "trajectory", "trajectory_manager", "vectors", "zobjects"
    ],
    package_data={
        "x11view": ["*.py"],
    },
    install_requires=[
        "numpy>=1.19.0",
        "python-xlib>=0.31",
    ],
    entry_points={
        "console_scripts": [
            "quickortep=quickORTEP:main",
        ],
    },
)