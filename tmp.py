class FileManager:
    def __init__(self, remote_name, local_master_dir, cloud_master_dir):
        self.remote_name = remote_name
        self.local_master_dir = local_master_dir
        self.cloud_master_dir = cloud_master_dir

    def convertCloudtoLocal(self, cloud_filename):
        relative_path = cloud_filename.replace(self.remote_name + self.cloud_master_dir, self.local_master_dir)
        return relative_path

    def convertLocaltoCloud(self, local_filename):
        relative_path = local_filename.replace(self.remote_name + self.local_master_dir, self.cloud_master_dir)
        return relative_path

    def uploadData (self, filename):
        cloud_path = self.convertLocaltoCloud(filename)
        return cloud_path

    def downloadData (self, filename):
        local_path = self.convertCloudtoLocal(filename)
        return local_path



my_fm = FileManager('remote:/', 'cloud_master/temp', 'local_master/temp')
