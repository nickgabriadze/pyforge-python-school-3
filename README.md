# Deployment to EC2 Instance



### Steps Taken

1. **Created `ssh-action.yml`**

   Created a GitHub Actions workflow file named `ssh-action.yml` to manage the deployment process via SSH.

2. **Started EC2 Instance**

   Logged into the AWS Management Console and manually started an EC2 instance to serve as the deployment target.

3. **Created Directory on EC2**

   SSH’d into the EC2 instance and created a directory named `projects` where the repository would be stored.

   ```bash
   mkdir -p /home/projects
4. **Cloned Repository into `projects` Directory**

   Cloned the repository containing the application code into the `projects` directory on the EC2 instance.
5. **Configured EC2 Environment**

   Configured the EC2 environment for the repository by setting up GitHub Secrets:
   - **`SSH_HOST`**: The public IP address of the EC2 instance.
   - **`SSH_KEY`**: The private key of the key-pair used for SSH access.

6. **Added Deployment Script**

   Added a script to the GitHub Actions workflow to automatically pull all the changes from the repository.

7. **Verified Deployment**

   After running the action, it successfully SSH’d into the EC2 instance and pulled the latest code from the repository, thus updating the project contents.
