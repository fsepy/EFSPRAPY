import csv
from typing import List


def csv_to_list_of_dicts(fp: str) -> List[dict]:
    """This is used to read csv file containing measurement data exported from Bluebeam.

    :param fp:
    :return:
    """
    with open(fp, newline='', encoding='utf-8-sig') as csvfile:
        # Read the CSV file
        reader = csv.reader(csvfile)

        # Get the headers (first row)
        headers = next(reader)

        # Initialize the result list
        result = []

        # Iterate through the CSV file and create the list of dictionaries
        for row in reader:
            row_dict = {header: value for header, value in zip(headers, row)}
            result.append(row_dict)

    return result


def csv_to_dict_of_lists(file_path):
    """
    Convert a CSV file to a dictionary of lists.

    Each column header in the CSV file becomes a key in the dictionary,
    and the items in that column become the values in the list.

    Args:
    file_path (str): The path to the CSV file.

    Returns:
    dict: A dictionary where keys are column headers and values are lists of column data.
    """
    with open(file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        result = {}

        # Initialize lists for each header
        for header in reader.fieldnames:
            result[header] = []

        # Populate the lists with column data
        for row in reader:
            for header in reader.fieldnames:
                result[header].append(row[header])

        return result


def dict_of_ndarray_to_csv(filename, data_dict, decimal_places=5):
    """
    Saves a dictionary of lists or numpy arrays to a CSV file,
    with optional numerical formatting for floats.

    Args:
        data_dict (dict): A dictionary where keys are column headers
                          and values are lists or numpy arrays of the same length.
        filename (str): The name of the CSV file to save.
        decimal_places (int, optional): If provided, float values will be
                                        formatted to this number of decimal places.
                                        Integers and other types are not affected
                                        by this formatting. Defaults to None.
    """
    if not data_dict:
        print("Warning: The dictionary is empty. An empty CSV file will be created.")
        with open(filename, 'w', newline='') as csvfile:
            pass  # Create an empty file
        return

    # Get the headers (keys of the dictionary)
    headers = list(data_dict.keys())  # Use list() to ensure a consistent order if needed

    # Get the number of rows. Assume all lists/arrays have the same length.
    # We can take the length of the first value.
    first_key = headers[0]
    num_rows = len(data_dict[first_key])

    # Optional: Add a check to ensure all lists/arrays have the same length
    for key, value in data_dict.items():
        if len(value) != num_rows:
            raise ValueError(f"Length mismatch: Value for key '{key}' has length {len(value)}, expected {num_rows}")

    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Write the header row
        writer.writerow(headers)

        # Write the data rows
        for i in range(num_rows):
            row_data = []  # List to hold the data for the current row
            for header in headers:
                value = data_dict[header][i]

                # Check if the value is a float and decimal formatting is requested
                if isinstance(value, float) and decimal_places is not None:
                    try:
                        # Format the float to the specified number of decimal places
                        # f-string formatting: :.nf where n is decimal_places
                        formatted_value = f"{value:.{decimal_places}f}"
                        row_data.append(formatted_value)
                    except ValueError:
                        # In case formatting fails for some unexpected reason
                        # Append the original value and let csv.writer handle it
                        row_data.append(value)
                # You could add an elif here if you want to specifically format integers
                # elif isinstance(value, int):
                #     # Example: format integers with leading zeros or specific width
                #     row_data.append(f"{value:05d}")
                else:
                    # For non-floats, integers (if not formatted above), None, strings, etc.
                    # Append the original value. csv.writer will convert it to its
                    # default string representation (e.g., 10 becomes "10", "hello" becomes "hello").
                    row_data.append(value)

            writer.writerow(row_data)
