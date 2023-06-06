## Deploying the site
pkgdown:::build_site_external()
setwd("docs")
rsconnect::deployApp(".",                       # the directory containing the content
                     appFiles = list.files(".", recursive = TRUE), # the list of files to include as dependencies (all of them)
                     appPrimaryDoc = "index.html",                 # the primary file
                     appName = "ECxHelper",                          # name of the endpoint (unique to your account on Connect)
                     appTitle = "Dose-Response Helper Functions",                         # display name for the content
                     account = "zhenglei_gao",                # your Connect username
                     server = "rsconnect.cs.sats.cloud"                    # the Connect server, see rsconnect::accounts()
)
