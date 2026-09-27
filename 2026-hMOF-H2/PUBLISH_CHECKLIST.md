# Before sending the repository link

1. Create an empty private GitHub repository named, for example, `hMOF-H2-reproducibility` (do not add a second README or license during creation).
2. From this directory, configure your Git identity if needed, then commit and push:

   ```bash
   git config user.name "Your Name"
   git config user.email "your.email@institution.edu"
   git commit -m "Prepare hMOF-H2 reproducibility package"
   git branch -M main
   git remote add origin https://github.com/USERNAME/hMOF-H2-reproducibility.git
   git push -u origin main
   ```

3. Upload the raw RASPA archives to approved shared storage. Replace `PASTE_SHARED_RAW_DATA_LINK_HERE` in `raw_manifest/README.md` with the real link and commit/push that change.
4. Grant the professor access to both the GitHub repository and the raw-data folder, or make both links accessible to the intended recipients.
5. Open the GitHub README in a private/incognito browser window and test that the repository link, raw-data link, and the validation commands are understandable before sending the email.

## Adding this project to the group repository

The group repository is organized as a multi-project archive. If that repository is the intended destination, copy this directory (without its `.git` directory) into one top-level project folder, for example `<publication-year>-hMOF-H2`. Work on a branch and submit the folder as one reviewable change; do not mix it with the unrelated historical directories. Run the commands in `README.md` from inside the new project folder.

Keep the multi-gigabyte RASPA archives in approved external storage and put the download link in `raw_manifest/README.md`. The GitHub folder should contain the curated summary data, normalized structures, the three metal-specific GCMC input bundles, screening/validation manifests, scripts, and checksums.

Do not upload the full local `hydrogen` working archive, the sibling archived `reproducibility_v5` directory, cluster credentials, private keys, or site-specific account files.
