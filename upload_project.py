from mobie.metadata import add_remote_project_metadata


# TODO
def add_metadata():
    bucket_name = "plankton-fibsem"
    endpoint = "https://s3.embl.de"
    add_remote_project_metadata('./data', bucket_name, endpoint)


def upload():
    pass


if __name__ == '__main__':
    add_metadata()
    # upload()
