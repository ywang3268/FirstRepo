import subprocess

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

    def uploadData (self, local_filename):
        cloud_path = self.convertLocaltoCloud(local_filename)
        subprocess.run(['rclone', 'copyto', local_filename, cloud_path])
        return cloud_path

    def downloadData (self, cloud_filename):
        local_path = self.convertCloudtoLocal(cloud_filename)
        subprocess.run(['rclone', 'copyto', cloud_filename, local_path])
        return local_path




my_fm = FileManager('remote:/', 'local_master/tmp', 'cloud_master/tmp')

