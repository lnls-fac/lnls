"""Repeating timer."""

import threading as _threading
import time as _time
import signal


class TimerError(Exception):
    pass


class Timer(object):
    """Repeating timer that calls function with the provided arguments every
    time interval is elapsed.
    """

    def __init__(self, interval, function, signal_handler=None, args=(), kwargs={}):
        """interval -- interval between function calls, in seconds
        function -- function to be called
        args -- tuple of arguments to be passed to function
        kwargs -- dictionary of arguments to be passed to function
        """
        self._interval = interval
        self._function = function
        self._args = args
        self._kwargs = kwargs
        self._is_running = False
        self._signal_handler = signal_handler
        
        if signal_handler:
            signal.signal(signal.SIGINT, signal_handler)

    @property
    def interval(self):
        return self._interval

    @interval.setter
    def interval(self, value):
        """value -- new interval in seconds"""
        if self._is_running:
            raise TimerError('timer is running')
        else:
            self._interval = value

    @property
    def args(self):
        return self._args

    @args.setter
    def args(self, value):
        if self._is_running:
            raise TimerError('timer is running')
        else:
            self._args = value

    @property
    def kwargs(self):
        return self._kwargs

    @kwargs.setter
    def kwargs(self, value):
        if self._is_running:
            raise TimerError('timer is running')
        else:
            self._kwargs = value

    @property
    def is_running(self):
        return self._is_running

    def start(self):
        """Start timer; first function call occurs after interval is elapsed.
        """
        if self._is_running:
            raise TimerError('timer is already running')
        else:
            self._is_running = True
            self._timer()

    def stop(self):
        if not self._is_running:
            raise TimerError('timer is not running')
        else:
            self._worker.cancel()
            self._worker.join(None)
            self._is_running = False

    def _timer(self):
        self._worker = _Worker(self._interval, self._function, self._args,
                               self._kwargs)
        self._worker.start()


class _Worker(_threading.Thread):

    def __init__(self, interval, function, args=[], kwargs={}):
        super(_Worker, self).__init__()
        self._interval = interval
        self._function = function
        self._args = args
        self._kwargs = kwargs
        self._event = _threading.Event()

    def cancel(self):
        self._event.set()

    def run(self):
        self._next_call = _time.time()
        while not self._event.is_set():
            self._next_call += self._interval
            self._event.wait(self._next_call-_time.time())
            if not self._event.is_set():
                self._function(*self._args, **self._kwargs)
