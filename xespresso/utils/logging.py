import logging

def get_logger(name: str = "xespresso", level: int = logging.INFO) -> logging.Logger:
    """
    Returns a reusable logger instance with stream output.

    Args:
        name (str): Name of the logger.
        level (int): Logging level (e.g. logging.INFO, logging.DEBUG).

    Returns:
        logging.Logger: Configured logger instance.
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)

    if not logger.hasHandlers():
        handler = logging.StreamHandler()
        formatter = logging.Formatter("[%(levelname)s] %(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    return logger
