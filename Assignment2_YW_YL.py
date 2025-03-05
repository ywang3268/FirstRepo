import os
import subprocess

class FileManager:
    def __init__(self, local_master_dir, cloud_master_dir, rclone_remote):
        self.local_master_dir = local_master_dir
        self.cloud_master_dir = cloud_master_dir
        self.rclone_remote = rclone_remote

    def convertCloudToLocal(self, filename):
        local_path = os.path.join(self.local_master_dir, filename)
        return local_path

    def convertLocalToCloud(self, filename):
        local_path = os.path.join(self.local_master_dir, filename)
        cloud_path = local_path.replace(self.local_master_dir, self.rclone_remote + ":" + self.cloud_master_dir)
        return cloud_path

    def uploadData(self, filename):
        cloud_path = self.convertLocalToCloud(filename).strip(filename)
        local_path = self.convertCloudToLocal(filename)
        print(f"Local file: {filename}")
        print(f"Cloud path: {cloud_path}")
        try:
            subprocess.run([
                "rclone", "copy",local_path, cloud_path
            ], check=True)
            print(f"Uploaded {filename} to {self.rclone_remote}:{cloud_path}")
        except subprocess.CalledProcessError as e:
            print(f"Error uploading {filename}: {e}")

    def downloadData(self, filename):
        local_path = self.convertCloudToLocal(filename).strip(filename)
        cloud_path = self.convertLocalToCloud(filename)
        print(f"Cloud file: {filename}")
        print(f"Local path: {local_path}")
        try:
            subprocess.run([
                "rclone", "copy", cloud_path, local_path], 
                check=True
                )
            print(f"Downloaded {self.rclone_remote}:{filename} to {local_path}")
        except subprocess.CalledProcessError as e:
            print(f"Error downloading {filename}: {e}")

if __name__ == "__main__":
    my_fm = FileManager(
        local_master_dir='/Users/eva/Desktop',  
        cloud_master_dir='/Python', 
        rclone_remote='dropbox' 
    )



    #Upload the file to Dropbox
    my_fm.uploadData("Sample1.txt")


    # Download the file from Dropbox
    my_fm.downloadData("Sample1.txt")