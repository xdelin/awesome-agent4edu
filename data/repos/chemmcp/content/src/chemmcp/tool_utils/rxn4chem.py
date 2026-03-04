import logging
from time import sleep
import os
from typing import List, Tuple

from rxn4chemistry import RXN4ChemistryWrapper  # type: ignore

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPToolInitError, ChemMCPApiNotFoundError


logger = logging.getLogger(__name__)


class RXN4Chem(BaseTool):
    """Wrapper for RXN4Chem functionalities."""

    base_url: str = "https://rxn.res.ibm.com"
    sleep_time: int = 5

    rxn4chem_chemistry_wrapper = None

    def __init__(self, rxn4chem_api_key=None, init=True, interface='text'):
        """Init object."""
        super().__init__(init, interface=interface)

        rxn4chem_api_key = os.environ.get("RXN4CHEM_API_KEY")
        if rxn4chem_api_key is None:
            raise ChemMCPApiNotFoundError("The RXN4Chem API key is not set. Please set the RXN4CHEM_API_KEY environment variable.")

        self.rxn4chem_api_key = rxn4chem_api_key
        if RXN4Chem.rxn4chem_chemistry_wrapper is None:
            RXN4Chem.rxn4chem_chemistry_wrapper = RXN4ChemistryWrapper(
                api_key=self.rxn4chem_api_key, base_url=RXN4Chem.base_url
            )
            RXN4Chem.rxn4chem_chemistry_wrapper.create_project('ChemMCP')
        self.rxn4chem = RXN4Chem.rxn4chem_chemistry_wrapper
        if init:
            if self.rxn4chem.project_id is None:
                raise ChemMCPToolInitError("The RXN4Chem project ID cannot be initialized.")
        logger.debug("RXN4Chem project ID: %s" % self.rxn4chem.project_id)

    @staticmethod
    def retry(times: int, exceptions, sleep_time: int = 5):
        """
        Retry Decorator.

        Retries the wrapped function/method `times` times if the exceptions
        listed in ``exceptions`` are thrown
        :param times: The number of times to repeat the wrapped function/method
        :type times: Int
        :param Exceptions: Lists of exceptions that trigger a retry attempt
        :type Exceptions: Tuple of Exceptions
        """

        def decorator(func):
            def newfn(*args, **kwargs):
                attempt = 0
                while attempt < times:
                    try:
                        sleep(sleep_time)
                        return func(*args, **kwargs)
                    except exceptions:
                        print(
                            "Exception thrown when attempting to run %s, "
                            "attempt %d of %d" % (func, attempt, times)
                        )
                        attempt += 1
                return func(*args, **kwargs)

            return newfn

        return decorator
