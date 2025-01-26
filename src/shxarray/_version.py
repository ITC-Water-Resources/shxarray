def version():
    try:
        import importlib
        return importlib.metadata.version('shxarray')
    except:
        return "unknown"
