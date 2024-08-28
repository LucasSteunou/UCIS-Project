class PathError(Exception):
    def __init__(self, message = 'No path found. Check if dep, end are in map, or if restricted areas allow a path.'):
        self.message = message
        super().__init__(self.message)

class WeightValueError(Exception):
    def __init__(self, message = 'Corridor Weight Value Error : Weight must be a number between 0 and 1.'):
        self.message = message
        super().__init__(self.message)
