class CustomSignal:
    def __init__(self):
        self._slots = []
    def connect(self, fn):
        self._slots.append(fn)
    def emit(self):
        for fn in self._slots:
            fn()
